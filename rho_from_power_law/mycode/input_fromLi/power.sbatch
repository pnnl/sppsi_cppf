#!/bin/bash
#SBATCH -A sppsi
### number of nodes
#SBATCH --nodes=3 --ntasks=64
### number of cores per node
####SBATCH --ntasks-per-node=24
#SBATCH -e power.e
#SBATCH -o power.o
### job name
#SBATCH -J power
#SBATCH -t 7-0:0:0
#SBATCH -p slurm

# Map each MPI process to core (-bycore).
mpirun -np 64 -bycore DAMASK_spectral --load power.load --geom input.geom
