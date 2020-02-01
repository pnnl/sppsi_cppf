from nice_plot import *

markers = ["b*", "rs", "gD"]

plt.figure()
for i in range(3):
    adir = "./power/postProc_rho/"
    forest = adir + "forest_inc%i.txt" % i
    forest = np.loadtxt(forest)
    
    # sum over 12 slip systems
    sforest = np.zeros(len(forest[0]))
    
    for j in range(len(forest)):
        sforest += forest[j]
        
    plt.plot(range(len(sforest)), sforest, markers[i], label = "forest t%i" % i)

# plt.yscale("log")
plt.legend(loc = 0)
plt.xlabel("Grid")
plt.ylabel("Dislocation density (m$^{-2}$)")
plt.savefig("rho_forest2.png")

plt.figure()
for i in range(3):
    adir = "./power/postProc_rho/"
    parallel = adir + "parallel_inc%i.txt" %i
    parallel = np.loadtxt(parallel)
    
    # sum over 12 slip systems
    sparallel = np.zeros(len(parallel[0]))
    
    for j in range(len(forest)):
        sparallel += parallel[j]
        
    plt.plot(range(len(sparallel)), sparallel, markers[i], label = "parallel t%i" %i)

# plt.yscale("log")
plt.legend(loc = 0)
plt.xlabel("Grid")
plt.ylabel("Dislocation density (m$^{-2}$)")
plt.savefig("rho_parallel2.png")