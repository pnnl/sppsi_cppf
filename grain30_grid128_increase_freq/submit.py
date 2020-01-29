import os
import numpy as np

# use 64 cores
# To update: JOBNAME, LOADFILE
sbatch_template = """#!/bin/bash
#SBATCH -A sppsi
### number of nodes
#SBATCH --nodes=3 --ntasks=64
### number of cores per node
####SBATCH --ntasks-per-node=24
#SBATCH -e JOBNAME.e
#SBATCH -o JOBNAME.o
### job name
#SBATCH -J JOBNAME
#SBATCH -t 3-0:0:0
#SBATCH -p slurm

# Map each MPI process to core (-bycore).
mpirun -np 64 -bycore DAMASK_spectral --load LOADFILE.load --geom 30grain_128grid.geom"""

# 30% gradient deformation
totalF = 0.3

# Load file of one slide cycle. TIME2 = 2 * TIME1. FDOT2 = -FDOT1
load = """fdot 0 FDOT1 0   0 0 0   0 0 *   stress * * *   * * *   * * -STRESS   time  TIME1  incs  100 freq 1
fdot 0 FDOT2 0  0 0 0  0 0 *  stress * * *  * * *  * * -STRESS time TIME2 incs 200"""

# strain rate from 1e-2 to 1e4 per second
fdot = np.array([1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3, 1e4])
# to label load file
strain_label = ['1en2', '1en1']
for i in range(2, len(fdot)):
	strain_label.append('%i' % fdot[i])

# time needed to achieve total strain
time = totalF / fdot

# 1st PK stress in MPa
stress = [20, 530]

for i in range(len(fdot)):
	for j in range(len(stress)):
		# back to parent dir from last submit
		os.system('cd /qfs/projects/sppsi/wkfu/grain30_grid128_increase_freq')
		# case name for label
		case = '30grain128grid_shearXY_fdot%s_p%i_power_law' % (strain_label[i], stress[j])
		os.system('mkdir %s' % case)
		
		# generate load file
		load_string = load.replace('FDOT1', '%.2e' % fdot[i])
		load_string = load_string.replace('TIME1', '%.2e' % time[i])
		load_string = load_string.replace('FDOT2', '-%.2e' % fdot[i])
		load_string = load_string.replace('TIME2', '%.2e' % (2.0 * time[i]))
		load_string = load_string.replace('STRESS', '%.2e' % (stress[j] * 1e6))
		with open('%s/%s.load' % (case, case), 'w') as f:
			f.write(load_string)
	
		# write sbatch file
		s = sbatch_template.replace('LOADFILE', case)
		s = s.replace('JOBNAME', 'f%sp%i' % (strain_label[i], stress[j]))
		with open('%s/%s.sbatch' % (case, case), 'w') as f:
			f.write(s)
	
		# copy unchanged files to running directory
		os.system('cp 30grain_128grid.geom material.config numerics.config Homogenization_None_Dummy.config Phase_Phenopowerlaw_Aluminum.config %s' % case)
		os.system('cd %s && sbatch %s.sbatch' % (case, case))
		
		
		
		
		
		
		
		
		



