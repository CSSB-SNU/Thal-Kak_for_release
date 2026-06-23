The Thal-Kak pipeline takes two YAML files: a data YAML and a model YAML.

## [Data YAML]
```
a3m:
- paired_path: str(AF3-like a3m Path) | null
  unpaired_path: str(AF3-like a3m Path) | null
  copy: int
  type: str(protein|dna|rna)
- ...

ligand (Optional):
- smiles: str(smiles) # exclusive with ccd
  copy: int
- ccd: str(CCD ID)    # exclusive with smiles
  copy: int
- ...

templates (Optional):
- path: str(cif Path)
  chain_template: str(Chain)
  chain_query:
  - str(Chain)
  - str(Chain)
  - ...

job_name: str
output_dir: str(Path)
seed:
- int
- int
- ...
```
### a3m field
- `paired_path`: path to the paired MSA for one sequence.
- `unpaired_path`: path to the unpaired MSA for one sequence.
- `copy`: how many copies of this sequence appear in the complex (e.g. `2` for an A2 complex).
- `type`: entity type of the sequence (e.g. `protein`).

For a query with two or more distinct sequences, prepare a separate paired
and unpaired MSA per sequence. Suppose an A1B1 stoichiometry gives the
following MSA:

```
Colabfold_MSA.a3m

#21,21  1,1   <-- ColabFold MSAs start with per-sequence lengths and copy counts; remove this line before use.
>101    102
AAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBB
|<-----A chain----->||<-----B chain----->|
>Uniref_something
AAAAAAAAAAAAAAAAAAAAABBBBBBBBBBBBBBBBBBBBB
|<-----A chain----->||<-----B chain----->|

>101
AAAAAAAAAAAAAAAAAAAAA---------------------
|<-----A chain----->||<-----B chain----->|
>Uniref_something
AAAAAAAAAAAAAAAAAAAAA---------------------
|<-----A chain----->||<-----B chain----->|

>102
---------------------BBBBBBBBBBBBBBBBBBBBB
|<-----A chain----->||<-----B chain----->|
>Uniref_something
---------------------BBBBBBBBBBBBBBBBBBBBB
|<-----A chain----->||<-----B chain----->|
```

For sequence A's paired MSA, drop the aligned B-chain columns and keep only
the A-chain part:
```
A_paired_MSA.a3m

>101    102
AAAAAAAAAAAAAAAAAAAAA
|<-----A chain----->|
>Uniref_something
AAAAAAAAAAAAAAAAAAAAA
|<-----A chain----->|
```
Likewise for sequence A's unpaired MSA, drop the aligned B-chain gaps:
```
A_unpaired_MSA.a3m

>101
AAAAAAAAAAAAAAAAAAAAA
|<-----A chain----->|
>Uniref_something
AAAAAAAAAAAAAAAAAAAAA
|<-----A chain----->|
```

For sequence B's paired MSA, drop the aligned A-chain columns:
```
B_paired_MSA.a3m

>101    102
BBBBBBBBBBBBBBBBBBBBB
|<-----B chain----->|
>Uniref_something
BBBBBBBBBBBBBBBBBBBBB
|<-----B chain----->|
```
And for sequence B's unpaired MSA, drop the aligned A-chain gaps:
```
B_unpaired_MSA.a3m

>102
BBBBBBBBBBBBBBBBBBBBB
|<-----B chain----->|
>Uniref_something
BBBBBBBBBBBBBBBBBBBBB
|<-----B chain----->|
```

### templates field
- `chain_template`: which chain of the template to use.
- `chain_query`: which query chain(s) the template is applied to.

Chain IDs are assigned in the order the sequences appear in the `a3m` field.
For example, if the first entity has copy 3 and the second has copy 2, the
first entity's copies become chains A, B, C and the second's become D, E.

## [Model YAML]
<details>
<summary><b>Boltz2</b></summary>
<div>

```
n_samples: int
subsample_msa: bool                     # default: False. True turns on boltz's native MSA subsampler (random 1024 rows per recycle)
constraints (Optional):
  - bond:
      atom1: [CHAIN_ID, RES_IDX, ATOM_NAME]
      atom2: [CHAIN_ID, RES_IDX, ATOM_NAME]
  - pocket:
      binder: CHAIN_ID
      contacts: [[CHAIN_ID, RES_IDX/ATOM_NAME], [CHAIN_ID, RES_IDX/ATOM_NAME]]
      max_distance: DIST_ANGSTROM		# default: 6.0
      force: false 						# default: false. If force is set to true, a potential will be used to enforce the pocket constraint
  - contact:
      token1: [CHAIN_ID, RES_IDX/ATOM_NAME]
      token2: [CHAIN_ID, RES_IDX/ATOM_NAME]
      max_distance: DIST_ANGSTROM		# default: 6.0
      force: false 						# default: false. If force is set to true, a potential will be used to enforce the contact constraint
```

Example:
```
n_samples: 5
subsample_msa: True
constraints:
  - pocket:
      binder: C
      contacts: [[A, 24], [A, 112]]
  - pocket:
      binder: D
      contacts: [[B, 24], [B, 112]]
```
See the [Boltz2 docs](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md) for details.

</div>
</details>

<details>
<summary><b>Chai1</b></summary>
<div>

```
num-trunk-samples:  		            # default: 1
num-diffn-samples:  		            # default: 5
num-diffn-timesteps:  		            # default: 200
recycle-msa-subsample:  		        # default: 0. 0 = full MSA; any positive value enables chai's per-recycle subsampler (kept-row count hardcoded to 4096 internally, so the magnitude is ignored)
num-trunk-recycles:  		            # default: 3
constraint-path: path to constraint     # default: null. If set, the constraints are used to enforce the structure prediction
use-esm-embeddings:  		            # default: True
fasta-names-as-cif-chains:  		    # default: False
template_hits_m8: path to m8 file       # default: null. If set, query/template alignment is read from this file instead of aligning manually when a template_path is provided in the data yaml
```

Example:
```
num-trunk-samples: 1
num-diffn-samples: 5
num-diffn-timesteps: 200
recycle-msa-subsample: 0
num-trunk-recycles: 3
constraint-path:
use-esm-embeddings: True
fasta-names-as-cif-chains: False
template_hits_m8:
```
See the [Chai1 docs](https://github.com/chaidiscovery/chai-lab/blob/main/README.md) for details.

</div>
</details>

Use the model YAML that matches the model you want to run. Start from the
example files and edit each field as needed.
