#!/bin/bash
#SBATCH --mail-user=lan.shi@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --job-name=perm
#SBATCH --time=4-23:59:00
#SBATCH --nodes=1    # comments allowed 
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --array=1-2400
#SBATCH --output=slurm_out/slurm-%A_%a.out

module load GCC/10.2.0 OpenMPI/4.0.5 R/4.0.5
R CMD BATCH --no-save --no-restore  alpha_permutation.R