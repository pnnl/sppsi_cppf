from read_damask import *

idir = "./30grain128grid_shearXY_fdot100_p20_power_law"
a = Dislocation(idir)
a.solver()