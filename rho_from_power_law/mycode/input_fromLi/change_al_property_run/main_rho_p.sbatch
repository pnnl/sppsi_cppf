#!/bin/bash
#SBATCH -A sppsi
### number of nodes
#SBATCH --nodes=1
### number of cores per node
#SBATCH --ntasks-per-node=1
#SBATCH -e rhop.e
#SBATCH -o rhop.o
### job name
#SBATCH -J rhop
#SBATCH -t 60:0:0
#SBATCH -p slurm


ipython main_rho.py ./pressure_shear
