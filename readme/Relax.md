# Relax

Relaxation is the final stage of the Thal-Kak pipeline. Given a directory of
decoy PDBs, it writes relaxed copies into `<decoy_dir>/relaxed/<method>/`. In
`full` pipeline mode this stage relaxes the per-job top-5 (selected by the
model's own confidence) in place.

## Where this stage sits

```
Structure  ──►  top-5 (top5/<job>/model_{1..5}.pdb)  ──►  [Relax]  ──►  top5/<job>/relaxed/<relax>/
```

## Methods (`--relax`)

| Option | Behavior |
|--------|----------|
| `none` | Pass-through. Copies each `*.pdb` to `<decoy_dir>/relaxed/none/<name>_unrelaxed.pdb`. |
| `amber` | OpenMM with the `amber19-all` force field + implicit solvent (OBC2 default; GBn2 or vacuum selectable in the script). Writes `<name>_relaxed_amber_<solvent>.pdb` (e.g. `_relaxed_amber_obc2.pdb`) plus per-decoy `<out_prefix>.energy.yaml` containing `E_init` / `E_min` / `E_final`, merged into `energies.yaml`. |

Both methods skip decoys whose output already exists, so re-running over a
partially-relaxed directory is cheap.

## CLI

```
thalkak relax --decoy_dir <dir of decoy PDBs> --relax amber
```

## Inputs / outputs

- **Input**: `<decoy_dir>/*.pdb`
- **Output**: `<decoy_dir>/relaxed/<method>/`
  - `<name>_relaxed_amber_<solvent>.pdb` (amber) or `<name>_unrelaxed.pdb` (none)
  - `energies.yaml`: `{ <decoy_basename>: { E_init, E_min, E_final, tool } }` (amber writes per-decoy `*.energy.yaml`, merged into `energies.yaml` by the orchestrator)

## What `amber` does

The `amber` method is adapted from AlphaFold2's relaxation (`alphafold.relax`),
with implicit solvent (instead of AF2's vacuum) and pLDDT-dependent coordinate
restraints. High-confidence
regions are held tight, and low-confidence regions move more. It runs in the
unified `thalkak` conda env (OpenMM is bundled there).

## Caveats

- Point at the flattened `common/` directory (or a top-5 directory), **not**
  the raw `seed_*/` subdirectories produced by some predictors.
- `relax=none` still writes files (with an `_unrelaxed` suffix) so downstream
  consumers see a consistent input contract.
