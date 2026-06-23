# MSA

MSA and template generation — stage 1 of the Thal-Kak pipeline. Takes a CASP-style FASTA and stoichiometry string, and produces the per-chain multiple sequence alignments, RNA/DNA chain handling, AF2 template hits, and the **data yaml** that every downstream stage consumes.

## Where this stage sits

```
FASTA + stoi  ──►  [MSA]  ──►  data yaml  ──►  Structure  ──►  Relax
```

## Methods (`--msa`)

| Option | Backend | What it does |
|--------|---------|--------------|
| `colab` | ColabFold (`colabfold_search`) via the configured env | MMseqs2 search against ColabFold DBs, then AF2 template lookup; the combined a3m is split per chain into paired / unpaired files |

RNA / DNA handling is automatic, regardless of `--msa`:
- Each FASTA record is checked against the nucleic-acid alphabet (`{A, C, G, T, U}`); records that match are routed out of the protein path.
- **RNA chains** → NHMMER-based MSA search against `rna_msa_db_dir`, output as a3m.
- **DNA chains** → no MSA; the FASTA is referenced directly in the data yaml.

## CLI

```
thalkak msa --msa colab \
            --seq examples/sample/T1201.fa \
            --stoi A1 \
            [--output_dir DIR]
```

## Inputs

- `--seq`: CASP FASTA. One record per distinct sequence, in chain order.
- `--stoi`: stoichiometry string, e.g. `A1`, `A2B1`, `A1B1C2`. `An` (literal `n`) marks an unknown copy count for chain `A` and is treated as `1`.
- `--output_dir`: optional override; defaults to the directory containing the FASTA.

## Outputs

Under `<output_dir>/msa/<msa_method>/`:

| File | Produced for |
|------|--------------|
| `<target>_paired_msa_chains_<x>.a3m` | each protein chain |
| `<target>_unpaired_msa_chains_<x>.a3m` | each protein chain |
| `<target>_na_<i>.a3m` | each RNA chain |
| `<target>_na_<i>.fa` | each DNA chain |
| `<pdb>_<chain>.cif` | AF2 template hits, when found |
| `method_log.yaml` | `{msa, seq, stoi, templates}` — inherited by Structure |

And one **data yaml** at `<output_dir>/<target>.yaml` matching the [data yaml schema](Structure.md#data-yaml-schema). The header reminds you to fill in `job_name`, `output_dir`, and `seed` before running structure prediction; in `thalkak full` mode these are filled in automatically.

## Skip-if-already-done

Each call writes / checks `method_log.yaml`. If `(msa, seq, stoi)` matches a previous run and `*.a3m` files already exist, MSA generation is skipped and the cached files are reused.

## Caveats

- Stoichiometry letters do not have to start at `A`; they are paired positionally with FASTA records.
- Chain IDs in downstream structure outputs are assigned in the order chains appear in the data yaml (see the templates-field section of [Structure.md](Structure.md)).
- The ColabFold raw a3m contains a header line like `#L1,L2  C1,C2` that must be removed when hand-preparing a3m. The auto path handles this for you.