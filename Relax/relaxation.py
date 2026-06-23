import argparse, subprocess, yaml, os, glob, shutil

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _merge_per_job_energies(relax_dir):
    """Aggregate per-job ``*.energy.yaml`` files (written by each relax worker)
    into a single ``energies.yaml`` keyed by the output basename. Done in the
    orchestrator after the parallel batch so workers don't race on a shared
    file."""
    energies_path = os.path.join(relax_dir, "energies.yaml")
    merged = {}
    if os.path.exists(energies_path):
        with open(energies_path) as f:
            merged = yaml.safe_load(f) or {}
    per_job_files = glob.glob(os.path.join(relax_dir, "*.energy.yaml"))
    for per_job in per_job_files:
        with open(per_job) as f:
            merged.update(yaml.safe_load(f) or {})
    if merged:
        with open(energies_path, "w") as f:
            yaml.dump(merged, f, sort_keys=False)
    for per_job in per_job_files:
        os.remove(per_job)


def relaxation(decoy_dir, relax_method):
    relax_dir = os.path.join(decoy_dir, "relaxed", relax_method)
    os.makedirs(relax_dir, exist_ok=True)

    match relax_method:
        case "none":
            print("No relaxation.")
            print(
                f"Copying unrelaxed structures from {decoy_dir} to relax directory..."
            )
            decoys = glob.glob(os.path.join(decoy_dir, "*.pdb"))
            for pdb in decoys:
                out_prefix = os.path.join(
                    relax_dir,
                    os.path.splitext(os.path.basename(pdb))[0] + "_unrelaxed",
                )
                if os.path.exists(f"{out_prefix}.pdb"):
                    print(
                        f"  Unrelaxed structure already exists: {out_prefix}.pdb, skipping."
                    )
                    continue
                shutil.copy(pdb, f"{out_prefix}.pdb")
            print("Copying complete.")

        case "amber":
            print("Running Amber (OpenMM, amber19-all + OBC2) relaxation...")
            decoys = glob.glob(os.path.join(decoy_dir, "*.pdb"))
            print(f"  found {len(decoys)} pdb files in {decoy_dir}")
            for pdb in decoys:
                out_prefix = os.path.join(
                    relax_dir,
                    os.path.splitext(os.path.basename(pdb))[0] + "_relaxed_amber_obc2",
                )
                if os.path.exists(f"{out_prefix}.pdb"):
                    print(
                        f"  Relaxed structure already exists: {out_prefix}.pdb, skipping."
                    )
                    continue
                subprocess.run(
                    f"python {ROOT}/Relax/script/amber/relax_amber.py "
                    f"-pdb_fn {pdb} -out_prefix {out_prefix} -implicit_solvent obc2",
                    shell=True,
                    check=True,
                )
            _merge_per_job_energies(relax_dir)
            print("Relaxation complete.")

    with open(os.path.join(relax_dir, "method_log.yaml"), "w") as f:
        yaml.dump({"relax": relax_method}, f, sort_keys=False)

    return relax_dir


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--decoy_dir",
        type=str,
        required=True,
        help="Path to directory containing decoy PDB files",
    )
    parser.add_argument(
        "--relax",
        type=str,
        required=True,
        choices=["none", "amber"],
        help="Relaxation method to use",
    )
    args = parser.parse_args()
    relaxation(args.decoy_dir, args.relax)
