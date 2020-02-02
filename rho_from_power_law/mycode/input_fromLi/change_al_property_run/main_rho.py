from read_damask import *

idir = sys.argv[1]
a = Dislocation(idir)
a.solver()