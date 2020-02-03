#!/bin/bash
#SBATCH -A sppsi
### number of nodes
#SBATCH --nodes=1
### number of cores per node
#SBATCH --ntasks-per-node=1
#SBATCH -e f1e4rho.e
#SBATCH -o f1e4rho.o
### job name
#SBATCH -J f1e4rho
#SBATCH -t 60:0:0
#SBATCH -p slurm


ipython main_rho_f1e4.py
