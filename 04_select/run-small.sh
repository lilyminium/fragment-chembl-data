#!/usr/bin/env bash
#SBATCH -J select
#SBATCH --array=0-36
#SBATCH -p standard
#SBATCH -t 16:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --account dmobley_lab
#SBATCH --output slurm-%x.%A-%a.out

# ===================== conda environment =====================
. ~/.bashrc
conda activate openff-nagl-test

date

export OE_LICENSE="/data/homezvol3/lilyw7/oe_license.txt"

TORSION_ID_ARRAY=(
    't126'
    't128'
    't8'
    't114'
    't112'
    't164'
    't18b'
    't81'
    't141b'
    't138a'
    't141a'
    't113'
    't167'
    't137'
    't7'
    't87a'
    't136'
    't102'
    't33'
    't141'
    't89'
    't103'
    't31a'
    't165'
    't141c'
    't55'
    't54'
    't73'
    't49'
    't42a'
    't101'
    't88'
    't158'
    't12'
    't154'
    't129'
    't30'
)

# Get the current index from the SLURM_ARRAY_TASK_ID environment variable
# and use it to select the corresponding TORSION_ID
TORSION_ID=${TORSION_ID_ARRAY[$SLURM_ARRAY_TASK_ID]}
echo "Running for TORSION_ID: $TORSION_ID"

python select-mols-with-low-torsions-small.py \
    -c "torsions-matching-${TORSION_ID}.csv"  \
    -s "selected-torsions-${TORSION_ID}.smi"  \
    -t $TORSION_ID


date
