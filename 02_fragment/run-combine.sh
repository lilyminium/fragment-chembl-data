#!/bin/bash
#SBATCH -J fragment-chembl
#SBATCH -p standard
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --account dmobley_lab
#SBATCH --export ALL
#SBATCH --constraint=fastscratch
#SBATCH --output slurm-%x.%A.out

source ~/.bashrc

conda activate yammbs


python combine-fragments-intermediate.py                 \
    -i fragments/chembl-35-fragments-single.smi     \
    -o output/chembl-35-fragments-combined.smi    \
    -np 64