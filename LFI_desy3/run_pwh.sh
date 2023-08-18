#!/bin/bash
#SBATCH -A des 
#SBATCH -C cpu 
#SBATCH -q regular 
#SBATCH -t 4:30:00 
#SBATCH --nodes=24
#SBATCH --ntasks=768


module load python
source activate /global/common/software/des/mgatti/py38_clone
cd /global/u2/m/mgatti/Mass_Mapping/peaks
srun python run_PWH.py




