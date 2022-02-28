#! /bin/bash
#SBATCH --job-name=ggm
#SBATCH --partition=long
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lorenzo.fabbri@isglobal.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=11
#SBATCH --mem=90gb
#SBATCH --output=/PROJECTES/HELIX/lorenzoF/logs/ggm.log

module purge > /dev/null 2>&1
module load lang/R

R CMD BATCH scripts_hpc.R
