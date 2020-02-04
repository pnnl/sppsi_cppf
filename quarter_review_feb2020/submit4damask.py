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
mpirun -np 64 -bycore DAMASK_spectral --load LOADFILE.load --geom GEOMFILE"""

tests = ["shear", "tensile"]
geom = [64, 128]
for test in tests:
	os.system("mkdir %s" % test)

# shear load file
load = """fdot * FDOT1 0   0 0 0   0 0 0   stress -STRESS * *   * * *   * * * time  TIME1  incs INCS freq 1"""

# strain rate from 1e-2 to 1e4 per second
fdot = np.array([1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3, 1e4])
# to label load file
strain_label = ['1en2', '1en1']
for i in range(2, len(fdot)):
	strain_label.append('%i' % fdot[i])

# 1st PK stress in MPa
stress = [20, 40]

# create root dirs
if "shear" not in os.listdir("./"):
	os.system("mkdir shear")
if "tensile" not in os.listdir("./"):
	os.system("mkdir tensile")

# # submit shear jobs
# for k in range(len(geom)):
# 	if geom[k] == 64:
# 				ffinal = 1.0
# 				incs = 400
# 				geom_file = "input_clean.geom"
# 				mat_file = "material64.config"
# 	else:
# 		ffinal = 0.3
# 		incs = 100
# 		geom_file = "30grain_128grid.geom"
# 		mat_file = "material128.config"
# 	
# 	time = ffinal / fdot
# 	
# 	for i in range(len(fdot)):
# 		for j in range(len(stress)):
# 			# make dir
# 			case = "s%i_f%s_p%i" %  (geom[k], strain_label[i], stress[j])
# 			folder = "./shear/%s" % case
# 			os.system("mkdir %s" % folder)
# 			
# 			# generate load file
# 			load_string = load.replace('FDOT1', '%.2e' % fdot[i])
# 			load_string = load_string.replace('STRESS', '%.2e' % (stress[j] * 1e6))
# 			load_string = load_string.replace('TIME1', '%.2e' % time[i])
# 			load_string = load_string.replace("INCS", "%i" % incs)
# 			with open('%s/%s.load' % (folder, case), 'w') as f:
# 				f.write(load_string)
# 				
# 			# write sbatch file
# 			s = sbatch_template.replace('LOADFILE', case)
# 			s = s.replace('JOBNAME', case)
# 			s = s.replace("GEOMFILE", geom_file)
# 			with open('%s/%s.sbatch' % (folder, case), 'w') as f:
# 				f.write(s)
# 			
# 			# back to parent dir from last submit
# # 			os.system('cd /qfs/projects/sppsi/wkfu/grain30_grid128_increase_freq')
# 
# 			# copy material file and change name
# 			os.system("cp %s %s/material.config" % (mat_file, folder))
# 			# copy other name-unchanged files
# 			unchanged_files = geom_file + " numerics.config Phase_Phenopowerlaw_Aluminum_mod.config"
# 			os.system('cp %s %s' % (unchanged_files, folder))
# 			os.system('cd %s && sbatch %s.sbatch' % (folder, case))
		
# submit tensile jobs
load = """fdot FDOT1 0 0   * * 0  * * *  stress * * *   0 0 * 0 0 0 time  TIME1  incs 100 freq 1"""
fdot = np.array([2e-3, 4e-3])
strain_label = ["2en3", "4en3"]
ffinal = 0.2
time = ffinal / fdot

for k in range(len(geom)):
	if geom[k] == 64:
				geom_file = "input_clean.geom"
				mat_file = "material64.config"
	else:
		geom_file = "30grain_128grid.geom"
		mat_file = "material128.config"

	for i in range(len(fdot)):
		# make dir
		case = "t%i_f%s" %  (geom[k], strain_label[i])
		folder = "./tensile/%s" % case
		os.system("mkdir %s" % folder)
		
		# generate load file
		load_string = load.replace('FDOT1', '%.2e' % fdot[i])
		load_string = load_string.replace('TIME1', '%.2e' % time[i])
		with open('%s/%s.load' % (folder, case), 'w') as f:
			f.write(load_string)
			
		# write sbatch file
		s = sbatch_template.replace('LOADFILE', case)
		s = s.replace('JOBNAME', case)
		s = s.replace("GEOMFILE", geom_file)
		with open('%s/%s.sbatch' % (folder, case), 'w') as f:
			f.write(s)
		
		# back to parent dir from last submit
# 			os.system('cd /qfs/projects/sppsi/wkfu/grain30_grid128_increase_freq')

		# copy material file and change name
		os.system("cp %s %s/material.config" % (mat_file, folder))
		# copy other name-unchanged files
		unchanged_files = geom_file + " numerics.config Phase_Phenopowerlaw_Aluminum_mod.config"
		os.system('cp %s %s' % (unchanged_files, folder))
		os.system('cd %s && sbatch %s.sbatch' % (folder, case))








