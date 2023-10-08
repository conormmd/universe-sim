#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=1G
#SBATCH --time=0:1:30

#SBATCH -p test.q

#SBATCH --mail-user=lgxn55@durham.ac.uk
#SBATCH --mail-type=ALL

module load python/3.6.8

python3 ./evolve_system.py