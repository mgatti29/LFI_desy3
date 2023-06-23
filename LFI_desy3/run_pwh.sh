#!/bin/bash 
#SBATCH -A des 
#SBATCH -C cpu 
#SBATCH -q regular 
#SBATCH -t 4:00:00 
#SBATCH --nodes=20
#SBATCH --ntasks=640

module load python
source activate py38
cd /global/u2/m/mgatti/Mass_Mapping/peaks
srun python run_PWH.py

