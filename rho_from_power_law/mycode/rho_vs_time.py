from nice_plot import *

idir = "./power/postProc_rho"

forest = []
parallel = []
mobile = []
gnd = []
ssd = []

nsteps = 10

for i in range(nsteps):
    ffile = idir + "/forest_inc%i.txt" % i
    pfile = idir + "/parallel_inc%i.txt" % i
    mfile = idir + "/mobile_inc%i.txt" % i
    gnd_file = idir + "/gnd_inc%i.txt" % i
    ssd_file = idir + "/ssd_inc%i.txt" % i
    
    f = np.loadtxt(ffile)
    p = np.loadtxt(pfile)
    m = np.loadtxt(mfile)
    g = np.loadtxt(gnd_file)
    s = np.loadtxt(ssd_file)
    
    forest.append(np.sum(f))
    parallel.append(np.sum(p))
    mobile.append(np.sum(m))
    gnd.append(np.sum(g))
    ssd.append(np.sum(s))
    
print("forest", forest)
print("gnd", gnd)
print("ssd", ssd)
print("mobile", mobile)
print("parallel", parallel)

x = range(len(forest))
plt.figure()
# plt.plot(x, forest, "b-", label = "forest")
# plt.plot(x, parallel, "g-", label = "parallel")
plt.plot(x, mobile, "r-", label = "mobile")
plt.plot(x, gnd, "k-", label = "GND")
plt.plot(x, ssd, "y-", label = "SSD")
# plt.yscale("log")
plt.legend(loc = 0)
plt.savefig("rho_vs_t.png")