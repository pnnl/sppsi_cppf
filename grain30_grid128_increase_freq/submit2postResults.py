import numpy as np
import os

"""Submit jobs to postProcess results, i.e., extract f,P from spectrcOut file"""

# one core for each case. To update JOBNAME
sbatch_template = """#!/bin/bash
#SBATCH -A sppsi
### number of nodes
#SBATCH --nodes=1
### number of cores per node
#SBATCH --ntasks-per-node=1
#SBATCH -e PRJOBNAME.e
#SBATCH -o PRJOBNAME.o
### job name
#SBATCH -J PRJOBNAME
#SBATCH -t 3-0:0:0
#SBATCH -p slurm

postResults --cr f,p *.spectralOut"""

# strain rate from 1e-2 to 1e4 per second
fdot = np.array([1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3, 1e4])
# to label load file
strain_label = ['1en2', '1en1']
for i in range(2, len(fdot)):
    strain_label.append('%i' % fdot[i])

# 1st PK stress in MPa
stress = [20, 530]

root_dir = '/qfs/projects/sppsi/wkfu/grain30_grid128_increase_freq'

for j in range(len(stress)):
    for i in range(len(fdot)):
        # case name for label
        case = '30grain128grid_shearXY_fdot%s_p%i_power_law' % (strain_label[i], stress[j])
        
        # cleanup
        os.system('cd %s/%s && rm -rf postProc' % (root_dir, case))
        
        # write sbatch file
        s = sbatch_template.replace('JOBNAME', 'f%sP%i' % (strain_label[i], stress[j]))
        with open('%s/%s/PR%s.sbatch' % (root_dir, case, case), 'w') as f:
            f.write(s)
        
        # submit job
        os.system('cd %s/%s && sbatch PR%s.sbatch' % (root_dir, case, case))