from read_damask_f1e4 import *

idir = "./30grain128grid_shearXY_fdot10000_p20_power_law"
a = Dislocation(idir)
a.solver()