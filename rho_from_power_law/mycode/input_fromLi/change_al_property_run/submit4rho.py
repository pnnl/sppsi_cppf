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

dirs = ["./pressure_shear", "./uniaxial_tensile2"]

for i in range(201):
    for j in range(2):
        idir = dirs[j]
        if j == 0:
            prefix = "p"
        else:
            prefix = "t"
        # process 100 steps
        s = sbatch_template.replace("PRJOBNAME", prefix + "i%i" % i)
        s = s.replace("DIR", idir)
        s = s.replace("INCS", "%i" % i)
        sbatch_fname = "rhoinc%i.sbatch" % i
        with open(idir + "/" + sbatch_fname, "w") as f:
            f.write(s)

        # submit
        os.system("cd %s && sbatch %s" % (idir, sbatch_fname))
