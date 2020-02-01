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

laws = ['isotropic', 'power', 'nonlocal']
prefile = "pressure_shear/postProc/input_power.txt"
tfile = "uniaxial_tensile2/postProc/input_power.txt"

pstrain = read_column(prefile, "2_ln(V)")
pstress = read_column(prefile, "2_Cauchy")
plt.figure()
plt.plot(pstrain, pstress * 1e-6, "k-")
plt.xlabel("Strain 12")
plt.ylabel("Cauchy stress 12 (MPa)")
plt.grid(1)
plt.savefig("ss_shear.png")

# fdot 1e-3 0 0  * * 0  * * *   stress * * *   0 0 *   0 0 0   time 200  incs 200 freq 1
pstrain = read_column(tfile, "1_ln(V)")
pstress = read_column(tfile, "1_Cauchy")
plt.figure()
plt.plot(pstrain, pstress * 1e-6, "k-")
plt.xlabel("Strain 11")
plt.ylabel("Cauchy stress 11 (MPa)")
plt.grid(1)
plt.savefig("ss_tensile.png")
