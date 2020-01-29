import numpy as np
import os

"""Submit jobs to calculate dislocation density"""

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
#SBATCH -t 3:0:0
#SBATCH -p short

ipython ../read_damask.py DIR INCS"""

root_dir = '/qfs/projects/sppsi/wkfu/grain30_grid128_increase_freq/30grain128grid_shearXY_fdot100_p20_power_law'


for i in range(101):
    # process 100 steps
    s = sbatch_template.replace("PRJOBNAME", "i%i" % i)
    s = s.replace("DIR", root_dir)
    s = s.replace("INCS", "%i" % i)
    sbatch_fname = "rhoinc%i.sbatch" % i
    with open(root_dir + "/" + sbatch_fname, "w") as f:
        f.write(s)
        
    # submit
    os.system("cd %s && sbatch %s" % (root_dir, sbatch_fname))