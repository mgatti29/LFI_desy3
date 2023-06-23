#!/bin/bash 
#SBATCH -A des 
#SBATCH -C cpu 
#SBATCH -q regular 
#SBATCH -t 4:00:00 
#SBATCH --nodes=20
#SBATCH --ntasks=1280

module load python
module load pytorch/1.10
source activate py38
pip install appdirs
cd /global/u2/m/mgatti/Mass_Mapping/peaks
srun python run_ST.py
