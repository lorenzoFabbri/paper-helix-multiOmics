#! /bin/bash
#SBATCH --job-name=ewaff
#SBATCH --partition=long
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=lorenzo.fabbri@isglobal.org
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=120gb
#SBATCH --output=/PROJECTES/HELIX/lorenzoF/logs/ewas.log

module purge > /dev/null 2>&1
module load lang/R

R CMD BATCH run_methylome2.R
