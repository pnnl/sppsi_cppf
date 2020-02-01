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

ipython ./main_process.py DIR INCS"""

root_dir = './30grain128grid_shearXY_fdot10000_p20_power_law'


for i in range(1, 101):
    # process 100 steps
    s = sbatch_template.replace("PRJOBNAME", "f1e4i%i" % i)
    s = s.replace("DIR", root_dir)
    s = s.replace("INCS", "%i" % i)
    sbatch_fname = "rhoinc%i.sbatch" % i
    with open(sbatch_fname, "w") as f:
        f.write(s)
        
    # submit
    os.system("sbatch %s" % sbatch_fname)