#!/bin/bash
#SBATCH --partition=GPU
#SBATCH --gpus=1
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-core=1
#SBATCH --mem=60000
#SBATCH --job-name="27K_SimpleImpute_Mean_1"
#SBATCH --array=0-4
#SBATCH --output=27K_Mean1-%a.txt
#SBATCH --error=27K_Mean1-%a.err
#SBATCH --mail-user=HENY4533@UWEC.EDU
#SBATCH --mail-type=END

module load python/3.9.2
module load python-libs/3.0

python methylation_dl_model.py $SLURM_ARRAY_TASK_ID