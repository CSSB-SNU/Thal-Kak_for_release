"""
AMBER relaxation via OpenMM with:
  - amber19-all force field (protein + NA + lipid + ions)
  - implicit solvent (OBC2 or GBn2)
  - flat-bottom harmonic restraint (pLDDT-dependent)
  - L-BFGS energy minimization
  - violation-informed iterative relax (AF2-style)

Restraint form (flat-bottom harmonic, isotropic 3D):
    d  = sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2)
    U  = 0.5 * k * (max(0, d - tol))^2

tol is the flat-bottom radius (Angstrom). Inside the well the atom is free;
outside it pays a quadratic penalty with spring constant k.
"""

import argparse
import io
import os
import sys
import time
import yaml
import numpy as np

import openmm
from openmm import unit
from openmm import app as openmm_app
from pdbfixer import PDBFixer


ENERGY = unit.kilocalories_per_mole
LENGTH = unit.angstroms

FF_PROTEIN = "amber19-all.xml"
IMPLICIT_SOLVENT_FF = {
    "obc2": "implicit/obc2.xml",
    "gbn2": "implicit/gbn2.xml",
    "none": None,
}


def get_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-pdb_fn", type=str, required=True)
    p.add_argument("-out_prefix", type=str, required=True)

    # restraint
    p.add_argument(
        "-lddt_cut",
        type=float,
        default=60.0,
        help="lDDT cutoff. Above -> hires params, below -> lowres params.",
    )
    p.add_argument(
        "-k_hires",
        type=float,
        default=10.0,
        help="spring constant for high-lDDT region (kcal/mol/A^2).",
    )
    p.add_argument(
        "-tol_hires",
        type=float,
        default=0.5,
        help="flat-bottom tolerance for high-lDDT region (A).",
    )
    p.add_argument(
        "-k_lowres",
        type=float,
        default=5.0,
        help="spring constant for low-lDDT region (kcal/mol/A^2).",
    )
    p.add_argument(
        "-tol_lowres",
        type=float,
        default=2.0,
        help="flat-bottom tolerance for low-lDDT region (A).",
    )
    p.add_argument(
        "-rst_set",
        type=str,
        default="non_hydrogen",
        choices=["non_hydrogen", "backbone", "c_alpha"],
    )

    # force field / solvent
    p.add_argument(
        "-implicit_solvent", type=str, default="obc2", choices=["obc2", "gbn2", "none"]
    )

    # minimizer
    p.add_argument("-max_iterations", type=int, default=2000)
    p.add_argument(
        "-tolerance",
        type=float,
        default=2.39,
        help="L-BFGS energy tolerance (kcal/mol).",
    )
    p.add_argument(
        "-max_outer_iterations",
        type=int,
        default=3,
        help="violation-informed iterations (AF2 default 20, ColabFold 3).",
    )
    p.add_argument(
        "-platform",
        type=str,
        default="CUDA",
        choices=["CUDA", "CPU", "OpenCL", "Reference"],
    )

    return p.parse_args()


def read_lddt(pdb_fn):
    """Read per-residue pLDDT from B-factor column (CA / C1' atoms)."""
    lddt = []
    with open(pdb_fn) as fp:
        for line in fp:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16].strip()
            if atom_name not in ("CA", "C1'"):
                continue
            lddt.append(float(line[60:66]))
    return lddt


def _strip_5prime_phosphate(topology, positions):
    """
    amber19-all expects 5'-OH terminus for RNA/DNA. If input PDB has a 5'
    phosphate (OP3 atom on the first residue of an NA chain), removing it
    is the simplest path to a templated system. Returns new (topology, positions).
    """
    NA_RES = {"A", "U", "G", "C", "T", "DA", "DT", "DG", "DC"}
    # To make a real 5'-OH terminus, the whole phosphate group must go:
    # P, OP1/O1P, OP2/O2P, OP3/O3P and any associated H atoms. Removing only
    # OP3 leaves a dangling P that no amber14 template matches.
    REMOVE_ATOM_NAMES = {
        "P",
        "OP1",
        "OP2",
        "OP3",
        "O1P",
        "O2P",
        "O3P",  # alternative naming
        "HOP2",
        "HOP3",
        "H5T",
    }

    # Identify atoms to drop
    drop_ids = set()
    for chain in topology.chains():
        first_residue = None
        for res in chain.residues():
            first_residue = res
            break
        if first_residue is None:
            continue
        if first_residue.name not in NA_RES:
            continue
        for atom in first_residue.atoms():
            if atom.name in REMOVE_ATOM_NAMES:
                drop_ids.add(atom.index)

    if not drop_ids:
        return topology, positions

    # Rebuild topology without those atoms
    new_top = openmm_app.Topology()
    new_pos_nm = []  # collect bare floats (nm) so we can wrap as one Quantity
    chain_map = {}
    res_map = {}
    atom_map = {}
    for chain in topology.chains():
        new_chain = new_top.addChain(chain.id)
        chain_map[chain] = new_chain
        for res in chain.residues():
            new_res = new_top.addResidue(res.name, new_chain, res.id, res.insertionCode)
            res_map[res] = new_res
            for atom in res.atoms():
                if atom.index in drop_ids:
                    continue
                new_atom = new_top.addAtom(atom.name, atom.element, new_res)
                atom_map[atom] = new_atom
                vec = positions[atom.index].value_in_unit(unit.nanometer)
                new_pos_nm.append([vec[0], vec[1], vec[2]])
    for bond in topology.bonds():
        a1, a2 = bond[0], bond[1]
        if a1 in atom_map and a2 in atom_map:
            new_top.addBond(atom_map[a1], atom_map[a2])

    new_pos = np.array(new_pos_nm, dtype=np.float64) * unit.nanometer

    print(f"  stripped {len(drop_ids)} 5'-phosphate atoms (OP3/HOP3/H5T)")
    return new_top, new_pos


def fix_pdb(pdb_fn):
    """Use PDBFixer to add missing atoms / hydrogens, return PDBFile."""
    fixer = PDBFixer(filename=pdb_fn)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)

    # Strip 5'-terminal phosphate BEFORE addMissingAtoms so PDBFixer matches the
    # 5'-OH template and adds H5T correctly. Doing this after would leave a
    # dangling P that no amber14 RNA template can match.
    new_top, new_pos = _strip_5prime_phosphate(fixer.topology, fixer.positions)
    fixer.topology = new_top
    fixer.positions = new_pos

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)

    buf = io.StringIO()
    openmm_app.PDBFile.writeFile(fixer.topology, fixer.positions, buf)
    buf.seek(0)
    return openmm_app.PDBFile(buf)


def select_restrained(atom, rst_set):
    """Return True if this atom should be restrained under given set."""
    if rst_set == "non_hydrogen":
        return atom.element.symbol != "H"
    elif rst_set == "backbone":
        return atom.name in (
            "N",
            "CA",
            "C",
            "O",
            "P",
            "OP1",
            "OP2",
            "O5'",
            "C5'",
            "C4'",
            "C3'",
            "O3'",
        )
    elif rst_set == "c_alpha":
        return atom.name in ("CA", "C1'")
    return False


def add_flat_bottom_restraint(
    system,
    ref_pdb,
    lddt_s,
    lddt_cut,
    k_hires,
    tol_hires,
    k_lowres,
    tol_lowres,
    rst_set,
    exclude_residues,
):
    """
    Flat-bottom harmonic restraint:
        d = sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2)
        U = 0.5 * k * (max(0, d - tol))^2

    k and tol are per-particle (set from per-residue lDDT).
    Returns number of restrained particles.
    """
    # +1e-8 nm^2 under the sqrt regularizes the gradient at d=0. Without it, an atom
    # sitting exactly on its restraint reference -- always true at iteration 0, since
    # the reference IS the (rounded) input structure -- gives d(d)/dx = (x-x0)/d =
    # 0/0 = NaN. step(d-tol)=0 there, but 0*NaN = NaN, which detonates minimizeEnergy
    # ("Particle coordinate is NaN"). The epsilon (d_min ~ 1e-4 A) is far below any
    # flat-bottom tolerance, so it is energetically negligible.
    energy_expr = (
        "0.5 * k * step(d - tol) * (d - tol)^2;"
        "d = sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2 + 1e-8)"
    )
    force = openmm.CustomExternalForce(energy_expr)
    force.addPerParticleParameter("k")
    force.addPerParticleParameter("tol")
    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    # OpenMM CustomExternalForce evaluates x, y, z in nanometers and energy
    # in kJ/mol. All per-particle parameters must be in those internal units.
    k_unit = unit.kilojoules_per_mole / (unit.nanometer**2)
    len_unit = unit.nanometer

    k_hires_u = (k_hires * ENERGY / (LENGTH**2)).value_in_unit(k_unit)
    k_lowres_u = (k_lowres * ENERGY / (LENGTH**2)).value_in_unit(k_unit)
    tol_hires_u = (tol_hires * LENGTH).value_in_unit(len_unit)
    tol_lowres_u = (tol_lowres * LENGTH).value_in_unit(len_unit)

    n_added = 0
    positions = ref_pdb.positions
    for atom in ref_pdb.topology.atoms():
        if atom.residue.index in exclude_residues:
            continue
        if not select_restrained(atom, rst_set):
            continue

        # per-residue lDDT (fall back to lowres if missing)
        res_idx = atom.residue.index
        if res_idx < len(lddt_s) and lddt_s[res_idx] > lddt_cut:
            k_val, tol_val = k_hires_u, tol_hires_u
        else:
            k_val, tol_val = k_lowres_u, tol_lowres_u

        # reference coordinates in nm (OpenMM internal unit)
        x0, y0, z0 = positions[atom.index].value_in_unit(len_unit)
        force.addParticle(atom.index, [k_val, tol_val, x0, y0, z0])
        n_added += 1

    system.addForce(force)
    return n_added


def build_system(pdb, implicit_solvent):
    """Build OpenMM System with amber19-all + optional implicit solvent."""
    ff_files = [FF_PROTEIN]
    if implicit_solvent != "none":
        ff_files.append(IMPLICIT_SOLVENT_FF[implicit_solvent])

    forcefield = openmm_app.ForceField(*ff_files)
    # For implicit-solvent systems, use CutoffNonPeriodic with 1 nm cutoff.
    # NoCutoff would force O(N^2) evaluations and OOM on systems > ~30k atoms.
    # For vacuum (implicit_solvent='none') keep NoCutoff to match AF2 behavior.
    if implicit_solvent == "none":
        nonbonded = openmm_app.NoCutoff
        kwargs = {}
    else:
        nonbonded = openmm_app.CutoffNonPeriodic
        kwargs = {"nonbondedCutoff": 1.0 * unit.nanometer}
    system = forcefield.createSystem(
        pdb.topology,
        nonbondedMethod=nonbonded,
        constraints=openmm_app.HBonds,
        rigidWater=True,
        **kwargs,
    )
    return system


from violations_torch import find_violations_torch


def overwrite_bfactor_with_lddt(pdb_str, lddt_s):
    """Rewrite B-factor with per-residue pLDDT. Keeps hydrogens."""
    lines_out = []
    prev_key = None
    i_res = -1
    for line in pdb_str.splitlines(keepends=True):
        if not line.startswith("ATOM"):
            lines_out.append(line)
            continue
        chain = line[21:22]
        resnum = line[22:26]
        icode = line[26:27]
        key = (chain, resnum, icode)
        if key != prev_key:
            i_res += 1
            prev_key = key
        bf = lddt_s[i_res] if i_res < len(lddt_s) else 0.0
        new_line = line[:60] + f"{bf:6.2f}" + line[66:]
        lines_out.append(new_line)
    return "".join(lines_out)


def run_minimize(pdb, ref_pdb, lddt_s, args, exclude_residues):
    """
    Single minimization pass.

    Args:
        pdb:     current PDBFile (topology + current positions, starting point)
        ref_pdb: original PDBFile (topology + AF2 positions, restraint reference)
                 must share atom order with `pdb`
    Returns (min_pdb_str, einit, efinal, next_pdb).
    """
    system = build_system(pdb, args.implicit_solvent)

    n_restrained = 0
    if args.k_hires > 0 or args.k_lowres > 0:
        n_restrained = add_flat_bottom_restraint(
            system,
            ref_pdb,
            lddt_s,
            args.lddt_cut,
            args.k_hires,
            args.tol_hires,
            args.k_lowres,
            args.tol_lowres,
            args.rst_set,
            exclude_residues,
        )

    integrator = openmm.LangevinIntegrator(
        0 * unit.kelvin, 1.0 / unit.picosecond, 1.0 * unit.femtosecond
    )
    try:
        platform = openmm.Platform.getPlatformByName(args.platform)
    except Exception as e:
        print(
            f"  WARNING: platform '{args.platform}' unavailable ({e}), "
            f"falling back to fastest available"
        )
        platform = None
    if platform is not None:
        sim = openmm_app.Simulation(pdb.topology, system, integrator, platform)
    else:
        sim = openmm_app.Simulation(pdb.topology, system, integrator)
    actual_platform = sim.context.getPlatform().getName()
    print(f"  platform: {actual_platform}")
    sim.context.setPositions(pdb.positions)

    state0 = sim.context.getState(getEnergy=True)
    einit = state0.getPotentialEnergy().value_in_unit(ENERGY)

    # OpenMM 8 minimizeEnergy expects force tolerance (kJ/mol/nm).
    # AF2 default 2.39 was historically passed as energy tolerance; in OpenMM 8
    # we pass it as kcal/mol/A and convert.
    force_tol = args.tolerance * ENERGY / LENGTH
    t_min0 = time.time()
    sim.minimizeEnergy(
        maxIterations=args.max_iterations,
        tolerance=force_tol,
    )
    print(f"  minimize: {time.time() - t_min0:.1f}s")

    state1 = sim.context.getState(getEnergy=True, getPositions=True)
    efinal = state1.getPotentialEnergy().value_in_unit(ENERGY)
    positions = state1.getPositions()

    buf = io.StringIO()
    openmm_app.PDBFile.writeFile(sim.topology, positions, buf)
    min_pdb_str = buf.getvalue()

    print(
        f"  restrained {n_restrained} atoms | "
        f"E_init={einit:.2f}  E_final={efinal:.2f} kcal/mol  "
        f"excluded_residues={len(exclude_residues)}"
    )

    # return reusable PDBFile for next iteration
    buf.seek(0)
    next_pdb = openmm_app.PDBFile(buf)
    return min_pdb_str, einit, efinal, next_pdb


def main():
    args = get_args()

    if os.path.exists(f"{args.out_prefix}.pdb"):
        print(f"Output exists: {args.out_prefix}.pdb -- skipping")
        return

    print(f"Input: {args.pdb_fn}")
    print(f"FF: {FF_PROTEIN} + implicit={args.implicit_solvent}")
    print(f"Restraint: flat-bottom harmonic, set={args.rst_set}")
    print(f"  hires (lDDT>{args.lddt_cut}): k={args.k_hires}, tol={args.tol_hires}")
    print(f"  lowres (lDDT<={args.lddt_cut}): k={args.k_lowres}, tol={args.tol_lowres}")

    t0 = time.time()

    # 1. fix & hydrogenate
    #    Note: amber19-all RNA templates assume 5'-OH terminus. If input has
    #    5'-phosphate (OP3), _strip_5prime_phosphate (called inside fix_pdb)
    #    removes it so the system can be templated. The minimized output will
    #    therefore represent a 5'-OH variant of the molecule.
    pdb = fix_pdb(args.pdb_fn)
    print(
        f"[fix_pdb] {time.time() - t0:.1f}s, "
        f"{pdb.topology.getNumAtoms()} atoms, "
        f"{pdb.topology.getNumResidues()} residues"
    )

    # Preserve original AF2 coordinates as the immutable restraint reference.
    # Re-parse the cleaned PDB so ref_pdb and the iterating `pdb` share atom order.
    ref_buf = io.StringIO()
    openmm_app.PDBFile.writeFile(pdb.topology, pdb.positions, ref_buf)
    ref_buf.seek(0)
    ref_pdb = openmm_app.PDBFile(ref_buf)

    # 2. read pLDDT from original input
    lddt_s = read_lddt(args.pdb_fn)
    print(
        f"[lddt] {len(lddt_s)} residues, "
        f"mean={np.mean(lddt_s):.1f}, "
        f"frac>{args.lddt_cut}={np.mean(np.array(lddt_s) > args.lddt_cut):.2f}"
    )

    # 3. violation-informed iterative minimization
    exclude_residues = set()
    einit_first = None
    for it in range(args.max_outer_iterations):
        print(f"\n[iter {it}]")
        min_pdb_str, einit, efinal, pdb = run_minimize(
            pdb, ref_pdb, lddt_s, args, exclude_residues
        )
        if einit_first is None:
            einit_first = einit

        t_v = time.time()
        # OpenMM CUDA and torch CUDA are independent installations. OpenMM may
        # be using its own bundled CUDA libs while torch was built against a
        # newer driver. Probe torch CUDA before committing; fall back to CPU
        # silently if it can't initialize. Violation detection is small enough
        # that CPU is fine.
        import torch as _torch

        device = "cpu"
        if _torch.cuda.is_available():
            try:
                _ = _torch.zeros(1, device="cuda")  # force lazy init
                device = "cuda"
            except Exception:
                device = "cpu"
        new_violations, vinfo = find_violations_torch(
            pdb,
            pdb.positions,
            tolerance_factor=12.0,
            clash_overlap_tolerance=1.5,
            device=device,
        )
        new_violations -= exclude_residues
        print(
            f"  violations: {len(new_violations)} res "
            f"(bond={vinfo['bond']}, between_clash={vinfo['between_clash']}, "
            f"within_clash={vinfo['within_clash']}) "
            f"[{time.time() - t_v:.1f}s]"
        )

        if not new_violations:
            print(f"  converged (no new violations)")
            break
        exclude_residues |= new_violations

    # 4. overwrite B-factor with pLDDT, drop hydrogens, then write output
    min_pdb_str = overwrite_bfactor_with_lddt(min_pdb_str, lddt_s)
    out_pdb = f"{args.out_prefix}.pdb"
    with open(out_pdb, "w") as f:
        f.write(min_pdb_str)

    # Energy reporting: E_init / E_min / E_final.
    # The AMBER pipeline has no separate relax step, so E_min == E_final.
    E_init = float(einit_first)
    E_min = float(efinal)
    E_final = float(efinal)
    tool_tag = f"amber_{args.implicit_solvent}"

    print(
        "ENERGY of INITIAL / MINIMIZED / RELAXED = %.3f / %.3f / %.3f (%s)"
        % (E_init, E_min, E_final, tool_tag)
    )

    # Per-job energy file (merged later by the orchestrator). Avoids
    # read-modify-write races when multiple workers run concurrently.
    e_file = f"{args.out_prefix}.energy.yaml"
    with open(e_file, "w") as fp:
        yaml.dump(
            {
                os.path.basename(args.out_prefix): {
                    "tool": tool_tag,
                    "E_init": float("%.6f" % E_init),
                    "E_min": float("%.6f" % E_min),
                    "E_final": float("%.6f" % E_final),
                }
            },
            fp,
            sort_keys=False,
        )

    print(f"[done] {time.time() - t0:.1f}s")
    print(f"  output:   {out_pdb}")
    print(f"  energies: {e_file}")
    print(f"  excluded {len(exclude_residues)} residues from restraint")


if __name__ == "__main__":
    main()

