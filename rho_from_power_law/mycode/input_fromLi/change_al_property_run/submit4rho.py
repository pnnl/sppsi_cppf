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

root_dir = ['./pressure_shear', "./uniaxial_tensile2"]

for j in range(len(root_dir)):
    idir = root_dir[j]
    if j == 0:
        prefix = "p"
    else:
        prefix = "t"
    for i in range(1, 201):
        # process 100 steps
        s = sbatch_template.replace("PRJOBNAME", prefix + "i%i" % i)
        s = s.replace("DIR", idir)
        s = s.replace("INCS", "%i" % i)
        sbatch_fname = prefix + "rhoi%i.sbatch" % i
        with open(sbatch_fname, "w") as f:
            f.write(s)
            
        # submit
        os.system("sbatch %s" % sbatch_fname)