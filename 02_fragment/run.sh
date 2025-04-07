#!/bin/bash
#SBATCH -J fragment-chembl
#SBATCH -p standard
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --constraint=fastscratch
#SBATCH --output slurm-%x.%A.out

source ~/.bashrc

conda activate yammbs


python fragment-on-rotatable-bonds.py               \
    -np 16   \
    -i ../01_download-initial/output/chembl_35.smi  \
    -o fragments/chembl-35-fragments.smi            \
    > logs/fragment-chembl-35.log 2>&1
