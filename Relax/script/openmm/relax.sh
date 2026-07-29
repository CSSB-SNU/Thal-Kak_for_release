#!/bin/bash
#SBATCH -p gpu
#SBATCH -J relax_openmm
#SBATCH --gres=gpu:A5000
#SBATCH --mem=16g
#SBATCH -c 4
#SBATCH -o log/%x_%j.log
#SBATCH -e log/%x_%j.err
#SBATCH --export=NONE

source ~/.bashrc
conda activate thalkak

# When run directly (bash/srun), cd to this script's own directory so the
# relative `python relax.py` resolves from any launch CWD. Under sbatch the
# script is copied to the spool dir (dirname $0 would point there), so leave CWD
# as SLURM_SUBMIT_DIR -- submit from this directory.
if [ -z "$SLURM_JOB_ID" ]; then
    cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")" || exit 1
fi

pdb_fn=$1
out_prefix=$2
lddt_cut=${3:-60.0}
implicit_solvent=${4:-obc2}   # obc2 | gbn2 | none
platform=${5:-CUDA}

# Normalize the input before relaxing (prochiral methyl names + C-terminal sp2
# carboxylate), matching the `thalkak relax` orchestrator, which validates up
# front. relax.py still keeps a post-PDBFixer OXT safety net for termini fixup.
here="$(cd "$(dirname "$(realpath "${BASH_SOURCE[0]}")")" && pwd)"
validated="${out_prefix}.validated_input.pdb"
python "$here/../validate.py" -pdb_fn "$pdb_fn" -out "$validated" -overwrite
pdb_fn="$validated"

python relax.py \
    -pdb_fn "$pdb_fn" \
    -out_prefix "$out_prefix" \
    -lddt_cut "$lddt_cut" \
    -implicit_solvent "$implicit_solvent" \
    -platform "$platform"

# example usage:
#   sbatch relax.sh examples/H1106/H1106.pdb examples/H1106/H1106_openmm 60 obc2
#   sbatch relax.sh examples/H1106/H1106.pdb examples/H1106/H1106_openmm 60 gbn2
