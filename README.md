<p align="center">
  <img src="assets/logo_light.png#gh-light-mode-only" width="600">
  <img src="assets/logo_dark.png#gh-dark-mode-only" width="600">
</p>

# Thal-Kak

A modular structure-prediction pipeline that runs several protein / nucleic-acid
structure predictors behind one interface. From an input FASTA it generates an
MSA, runs the chosen predictor across multiple seeds, selects the top-5 models
by the predictor's own confidence, and relaxes them.

## Pipeline overview

Three modular stages, each runnable on its own through the `thalkak` CLI:

- **MSA generation** — turns the input FASTA + stoichiometry into a normalized
  `data.yaml` (ColabFold for protein, NHMMER for nucleic acids).
- **Structure prediction** — takes `data.yaml` and a model config, dispatches to
  one of the supported backends, and flattens every output into a shared
  `common/` directory (decoy PDBs + PAE/pLDDT plots + a confidence CSV).
- **Relaxation** — relaxes the selected top-5 models and writes refined PDBs.

`thalkak full` composes MSA → Structure → top-5 selection → Relax.

### Available options
| Stage | Options |
|-------|---------|
| MSA (`--msa`) | `colab` |
| Structure (`--structure`) | `boltz2`, `chai1`, `protenix`, `esmfold2` |
| Relaxation (`--relax`) | `none`, `amber` |

Per-stage docs:
- [readme/MSA.md](readme/MSA.md) — ColabFold MSA + templates; RNA chains routed through NHMMER
- [readme/Structure.md](readme/Structure.md) — Boltz-2 / Chai-1 / Protenix / ESMFold2 runners, data/model yaml schemas
- [readme/Relax.md](readme/Relax.md) — OpenMM Amber relaxation with pLDDT-weighted restraints

### Data flow in `full` mode
```
FASTA + stoi
    │
    ▼
  [MSA]   ──►  msa/<method>/*.a3m  +  method_log.yaml  +  <target>.yaml (data yaml)
    │
    ▼
[Structure]  ──►  structure/<model>_results_<target>_<job>/common/
                       ├── *.pdb
                       ├── *_results_summary.csv   (confidence)
                       ├── *.png
                       └── method_log.yaml         (msa + structure)
    │
    ▼
[top-5]  ──►  top5/<job>/model_{1..5}.pdb + method_log.yaml
    │         (ranked by the model's own confidence; bond-length validated)
    ▼
  [Relax]  ──►  top5/<job>/relaxed/<relax>/   (+ energies.yaml)
```

`method_log.yaml` is threaded through every stage — each appends its own choice
(`msa`, `structure`, `relax`) so the output carries a provenance record.

RNA / DNA targets are auto-detected from the FASTA alphabet: RNA chains go
through an NHMMER-based MSA search, DNA chains pass through as FASTA.

## Usage

> Activate the environment first: `conda activate thalkak`

### Full pipeline
```
thalkak full --msa colab --structure boltz2 --relax amber \
    --seq examples/sample/T1201.fa --stoi A1
```

Additional `full` options:

| Flag | Default | Purpose |
|------|---------|---------|
| `--model_config <path>` | `examples/{structure}.yaml` | Override the per-model config yaml |
| `--n_seed <int>` | `5` | Number of seeds for structure prediction |
| `--seed_start <int>` | `1` | First seed value (seeds are `seed_start .. seed_start + n_seed - 1`) |
| `--base_dir <path>` | dir of `--seq` | Output root |
| `--ligand_yaml <path>` | none | YAML with a `ligand:` list merged into the data yaml. Per backend: `boltz2` accepts `ccd` + `smiles`; `chai1` accepts `smiles`; `protenix` accepts `ccd` + `smiles` + `file`. |

### Individual stages
```
# MSA only
thalkak msa --msa colab --seq examples/sample/T1201.fa --stoi A1

# Structure only (fill job_name / output_dir / seed in the data yaml first)
thalkak structure --model boltz2 \
    --data_config examples/sample/T1201.yaml \
    --model_config examples/boltz2.yaml

# Relax only
thalkak relax --decoy_dir <dir of decoy PDBs> --relax amber
```

## Install

```
bash install.sh
```

This creates the unified `thalkak` conda environment (Boltz-2, Chai-1,
Protenix, ESMFold2, ColabFold MSA, and OpenMM/Amber relaxation all in one env)
from `environment.yml` plus the `--no-deps` packages in
`requirements-nodeps.txt`. See the comments at the top of `install.sh` for the
exact steps.

### Configuration
Edit `config.yaml`:

| Key | Used by |
|-----|---------|
| `rna_msa_db_dir` | Local RNA MSA database (only needed for RNA targets). |

```yaml
rna_msa_db_dir: /path/to/rna_msa_db
```

Each structure-prediction model downloads its own weights to its default cache
location on first run (e.g. `~/.boltz` for Boltz-2, `~/.cache/huggingface` for
ESMFold2).

The vendored model sources live under `Structure/submodules/` and were pulled in
as git subtrees; see [readme/subtrees.yaml](readme/subtrees.yaml) for their
upstream origins.
