import os
from nice_plot import *

def read_column(afile, label):
	"""Return a column with label from the damask output afile."""
	with open(afile, 'r') as f:
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
	data = np.loadtxt(afile, skiprows = label_line + 1)
	return data[:, colID]

markers = ['b-', 'g--', 'r-.*', 'y:', 'm-', 'c--', 'k-.']

prefile = "pressure_shear/postProc/input_power.txt"
tfile = "uniaxial_tensile2/postProc/input_power.txt"

# fdot * 1e-3 0   1e-3 0 0   0 0 0   stress -2e7 * *   * * *   * * *   time 200  incs 200 freq 1
pstrain = read_column(prefile, "2_ln(V)")
pstress = read_column(prefile, "2_Cauchy")

# print final strain
print("Pressure load, final strain 12 = %.2f" % pstrain[-1])

# stored in txts are the sum of dislocation densities in 64^3 grids
# divided by 64**3 to compute averaged rho
ssd = np.loadtxt("pressure_shear/postProc_rho/ssd_sum.txt") / 64 ** 3
gnd = np.loadtxt("pressure_shear/postProc_rho/gnd_sum.txt") / 64 ** 3
plt.figure()
plt.plot(pstrain, pstress * 1e-6, "k-", label = "stress")
plt.xlabel("Strain 12")
plt.ylabel("Cauchy stress 12 (MPa)")
plt.legend(loc = 0)
plt.twinx()
plt.plot(pstrain[1 : ], ssd, "b--", label = "SSD")
plt.plot(pstrain[1 : ], gnd, "g--", label = "GND")
plt.ylabel("Dislo. density (m$^{-2}$)")
plt.yscale("log")
plt.legend(loc = "lower left", bbox_to_anchor = (0.5, .005))
plt.savefig("ss_shear.png")

# fdot 1e-3 0 0  * * 0  * * *   stress * * *   0 0 *   0 0 0   time 200  incs 200 freq 1
pstrain = read_column(tfile, "1_ln(V)")
pstress = read_column(tfile, "1_Cauchy")
ssd = np.loadtxt("uniaxial_tensile2/postProc_rho/ssd_sum.txt") / 64 ** 3
gnd = np.loadtxt("uniaxial_tensile2/postProc_rho/gnd_sum.txt") / 64 ** 3

# print final strain
print("Tensile load, final strain 11 = %.2f" % pstrain[-1])

plt.figure()
plt.plot(pstrain, pstress * 1e-6, "k-", label = "stress")
plt.xlabel("Strain 11")
plt.ylabel("Cauchy stress 11 (MPa)")
plt.legend(loc = 0)
plt.twinx()
plt.plot(pstrain[1 : ], ssd, "b--", label = "SSD")
plt.plot(pstrain[1 : ], gnd, "g--", label = "GND")
plt.ylabel("Dislo. density (m$^{-2}$)")
plt.yscale("log")
plt.legend(loc = "lower left", bbox_to_anchor = (0.5, .005))
plt.savefig("ss_tensile.png")
