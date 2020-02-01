from read_damask import *

idir = sys.argv[1]
inc = int(sys.argv[2])
a = Dislocation(idir)
a.post_process(inc)
