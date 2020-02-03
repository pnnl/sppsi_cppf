from nice_plot import *
import os

def read_column(txt, label):
        """Return a column with label from the damask output txt file."""
        with open(txt, 'r') as f:
            for i, line in enumerate(f):
                split_line = line.split()
                if split_line[0] == "inc":  # the label line
                    label_line = i
                    # the label row
                    for j in range(len(split_line)):
                        if split_line[j] == label:
                            colID = j
                            break
                    break
        data = np.loadtxt(txt, skiprows = label_line + 1, usecols = colID)
        return data
    
dirs = ["30grain128grid_shearXY_fdot100_p20_power_law", "30grain128grid_shearXY_fdot10000_p20_power_law"]

# strain rate from 1e-2 to 1e4 per second
fdot = np.array([1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3, 1e4])
# to label load file
strain_label = ['1en2', '1en1']
for i in range(2, len(fdot)):
    strain_label.append('%i' % fdot[i])

markers = ['b-', 'g-', 'r-', 'y-', 'm-', 'c-', 'k-']

# fdot shown in legend
fdot_legend = [r'$10^{-2}$', r'$10^{-1}$', '$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$']

# store strains for dislocation density plot
strains = []
# ss curves for different strain rates
plt.figure()
for i in range(len(fdot)):
    idir = '30grain128grid_shearXY_fdot%s_p20_power_law/postProc' % strain_label[i]
    txt = idir + "/" + os.listdir(idir)[0]
    strain = read_column(txt, "Mises(ln(V))")
    # used to plot dislocation density evolution
    strains.append(strain[1 : 101])
    stress = read_column(txt, "Mises(Cauchy)")
    if i == 0:
        print("length of strain is %i" % len(strain))
    plt.plot(strain[0 : 101], stress[0 : 101] * 1e-6, markers[i], label = fdot_legend[i])
plt.legend(loc = 0, ncol = 3)
plt.xlabel('Strain')
plt.ylabel('Stress (MPa)')
plt.text(0.12, 50, r'strain rate (s$^{-1}$)')
plt.ylim(ymax = 240)
plt.grid()
plt.savefig('ss_p20v2.png')

print("The final strain is %.2f" % strains[0][-1])

# dislocation density vs strain
plt.figure()
for i in range(len(fdot)):
    # we only have dislocation density of fdot = 100 and 10000
    if strain_label[i] in ["100"]:
        strain = strains[i]
        idir = '30grain128grid_shearXY_fdot%s_p20_power_law/postProc_rho' % strain_label[i]
        # stored in the txts are the sum of dislocation densities in the 128**3 grids, calculate its average over the entire volume
        ssd = np.loadtxt(idir + "/ssd_sum.txt") / 128 ** 3
        gnd = np.loadtxt(idir + "/gnd_sum.txt") / 128 ** 3
        plt.plot(strain, ssd, markers[i][0] + "-", label = r"SSD %s s$^{-1}$ strain rate" % fdot_legend[i])
        plt.plot(strain, gnd, markers[i][0] + ":", label = r"GND %s s$^{-1}$ strain rate" % fdot_legend[i])
plt.xlabel('Strain')
plt.ylabel("Dislo. density (m$^{-2}$)")
plt.grid()
plt.legend(loc = 0)
plt.yscale("log")
plt.savefig("dislo_vs_strain.png")
        
        
