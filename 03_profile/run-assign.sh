#!/usr/bin/env bash
#SBATCH -J batch
#SBATCH -p standard
#SBATCH -t 16:00:00
#SBATCH --nodes=1
#SBATCH --tasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=24gb
#SBATCH --account dmobley_lab
#SBATCH --output slurm-%x.%A.out

# ===================== conda environment =====================
. ~/.bashrc
conda activate openff-nagl-test

export OE_LICENSE="/data/homezvol3/lilyw7/oe_license.txt"

# python assign-torsions-and-profile.py                           \
#         --n-workers                     300                     \
#         --worker-type                   "slurm"                 \
#         --batch-size                    5                       \
#         --memory                        4                       \
#         --walltime                      480                     \
#         --queue                         "free"                  \
#         --conda-environment             "openff-nagl-test"      \
#     -i        "../02_fragment/test.smi"                         \
#     -o        "ff-parameters"


python assign-parameters-by-file.py                             \
        --n-workers                     300                     \
        --worker-type                   "slurm"                 \
        --batch-size                    1                       \
        --memory                        16                      \
        --walltime                      480                     \
        --queue                         "free"                  \
        --conda-environment             "openff-nagl-test"      \
    -i        "../02_fragment/intermediate"                     \
    -o        "ff-parameters"

