#!/bin/bash
#SBATCH -A sppsi
### number of nodes
#SBATCH --nodes=1
### number of cores per node
#SBATCH --ntasks-per-node=1
#SBATCH -e rhot.e
#SBATCH -o rhot.o
### job name
#SBATCH -J rhot
#SBATCH -t 60:0:0
#SBATCH -p slurm


ipython main_rho.py ./uniaxial_tensile2
