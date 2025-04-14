#!/usr/bin/env bash
#SBATCH -J select
#SBATCH --array=0-30
#SBATCH -p standard
#SBATCH -t 16:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128gb
#SBATCH --account dmobley_lab
#SBATCH --output slurm-%x.%A-%a.out

# ===================== conda environment =====================
. ~/.bashrc
conda activate openff-nagl-test

date

export OE_LICENSE="/data/homezvol3/lilyw7/oe_license.txt"

PARAMETER_ID_ARRAY=(
    # bonds < 50 representations
    'b49'
    'b78'
    'b55'
    'b82'
    'b80'
    'b40'
    'b29'
    'b81'
    'b63'
    'b44'
    'b76'
    # bonds where mm optimization shifts substantial points from qm
    'b13a'
    'b15'
    'b23'
    'b24'
    'b42'
    'b47'
    'b54'
    'b59'
    'b62'
    'b64'
    'b67'
    'b65'
    'b69'
    'b74'
    'b77'

    # angles < 50 representations
    'a36'
    'a16'
    'a7'

    # impropers with wide variation
    'i4'
    'i5'
)

# Get the current index from the SLURM_ARRAY_TASK_ID environment variable
# and use it to select the corresponding TORSION_ID
PARAMETER_ID=${PARAMETER_ID_ARRAY[$SLURM_ARRAY_TASK_ID]}
echo "Running for PARAMETER_ID: $PARAMETER_ID"

python select-mols-valence.py \
    -c "parameters-matching-${PARAMETER_ID}.csv"  \
    -s "selected-parameters-${PARAMETER_ID}.smi"  \
    -p $PARAMETER_ID


date
