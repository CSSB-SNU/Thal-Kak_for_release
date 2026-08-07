<p align="center">
  <img src="assets/logo_light.png#gh-light-mode-only" width="600">
  <img src="assets/logo_dark.png#gh-dark-mode-only" width="600">
</p>

# Thal-Kak

A modular structure-prediction pipeline that runs multiple protein and nucleic-acid
structure predictors through a single interface. Starting from an input FASTA, it
generates an MSA, runs the selected predictor across multiple seeds, selects the
top-5 models based on the predictor's own confidence score, and relaxes them.

**Try it in Colab (protein targets):**
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/CSSB-SNU/Thal-Kak_for_release/blob/main/Thalkak.ipynb)

> Click the badge, choose a **GPU** runtime (`Runtime → Change runtime type → GPU`),
> then `Runtime → Run all`. The environment installs in ~3 min (pixi); the
> first run also downloads model weights (a few minutes).

## Pipeline overview

<div align="center">
  <img alt="Thal-Kak pipeline scheme" src="assets/Thal-kak_scheme.png" width="800">
</div>

The pipeline consists of three modular stages, each of which can be run independently through the `thalkak` CLI:

- **MSA generation:** converts the input FASTA and stoichiometry into a normalized
  `data.yaml` file using ColabFold for proteins and NHMMER for nucleic acids.
- **Structure prediction:** takes `data.yaml` and a model config, dispatches to
  one of the supported backends, and flattens every output into a shared
  `common/` directory (decoy PDBs + PAE/pLDDT plots + a confidence CSV).
- **Relaxation:** relaxes the selected top-5 models and writes refined PDBs.

`thalkak full` composes MSA → Structure → top-5 selection → Relax.

### Available options
| Stage | Options |
|-------|---------|
| MSA (`--msa`) | `colab` |
| Structure (`--structure`) | `boltz2`, `chai1`, `protenix`, `esmfold2` |
| Relaxation (`--relax`) | `none`, `openmm` |

Per-stage documentation:
- [readme/MSA.md](readme/MSA.md): ColabFold MSA and templates; RNA chains are routed through NHMMER.
- [readme/Structure.md](readme/Structure.md): Boltz-2, Chai-1, Protenix, and ESMFold2 runners, plus data/model YAML schemas.
- [readme/Relax.md](readme/Relax.md): OpenMM all-atom relaxation with pLDDT-weighted restraints.

### Data flow in `full` mode
```
  FASTA + stoichiometry
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
    │         (ranked by the model's own confidence; bond length validated)
    ▼
 [Relax]  ──►  top5/<job>/relaxed/<relax>/   (+ energies.yaml)
```

`method_log.yaml` is threaded through every stage. Each stage appends its own
choice (`msa`, `structure`, `relax`), so the output carries a provenance record.

RNA and DNA targets are automatically detected from the FASTA alphabet: RNA chains go
through an NHMMER-based MSA search, while DNA chains pass through as FASTA.

## Install

```
bash install.sh
```

This creates the unified `thalkak` conda environment (Boltz-2, Chai-1,
Protenix, ESMFold2, ColabFold MSA, and OpenMM relaxation all in one env)
from `environment.yml` plus the `--no-deps` packages in
`requirements-nodeps.txt`. See the comments at the top of `install.sh` for the
exact steps.

### RNA MSA database (only for RNA targets)

RNA targets need a local RNA MSA database (Rfam + RNAcentral). Build it once
with the provided script, which downloads, clusters, and indexes the databases
into `MSA/script/RNA_MSA_search/db` — exactly where the pipeline looks for them:

```bash
conda activate thalkak
cd MSA/script/RNA_MSA_search
bash prepare_db.sh
```

The tools it uses (`mmseqs`, HMMER) are already in the conda environment. The
RNAcentral download and clustering is large and can take hours. Skip this step
entirely if you only run protein targets.

Each structure-prediction model downloads its own weights to its default cache
location on first run (e.g. `~/.boltz` for Boltz-2, `~/.cache/huggingface` for
ESMFold2).

Protenix is an exception: its official `protenix-v2` checkpoint endpoint
currently returns HTTP 403 for public requests. When the checkpoint is missing,
the pipeline downloads it from a community mirror and verifies it against a
pinned SHA-256 before use (a mismatch aborts the run). To supply your own copy
instead, place `protenix-v2.pt` in the protenix checkpoint directory
(`$PROTENIX_CHECKPOINT_DIR` if set, otherwise
`Structure/submodules/protenix/checkpoint/`).

The vendored model sources live under `Structure/submodules/` and were pulled in
as git subtrees; see [readme/subtrees.yaml](readme/subtrees.yaml) for their
upstream origins.

## Usage

> Activate the environment first: `conda activate thalkak`

### Full pipeline
```
thalkak full --msa colab --structure boltz2 --relax openmm \
    --seq examples/sample/T1201.fa --stoi A1
```

Additional `full` options:

| Flag | Default | Purpose |
|------|---------|---------|
| `--model_config <path>` | `examples/{structure}.yaml` | Override the per-model config yaml |
| `--n_seed <int>` | `5` | Number of seeds for structure prediction |
| `--seed_start <int>` | `1` | First seed value (seeds are `seed_start .. seed_start + n_seed - 1`) |
| `--base_dir <path>` | dir of `--seq` | Output root |
| `--ligand_yaml <path>` | none | YAML file with a `ligand:` list merged into the data yaml. Per backend: `boltz2` accepts `ccd` + `smiles`; `chai1` accepts `ccd` + `smiles`; `protenix` accepts `ccd` + `smiles` + `file`; `esmfold2` accepts `ccd` + `smiles`. |

### Individual stages
```
# MSA only
thalkak msa --msa colab --seq examples/sample/T1201.fa --stoi A1

# Structure only
#   The data yaml (examples/sample/T1201.yaml) is produced by the `msa` step
#   above (written next to the FASTA). Fill in job_name / output_dir / seed in
#   it before running structure prediction.
thalkak structure --model boltz2 \
    --data_config examples/sample/T1201.yaml \
    --model_config examples/boltz2.yaml

# Relax only
thalkak relax --decoy_dir <dir of decoy PDBs> --relax openmm
```
