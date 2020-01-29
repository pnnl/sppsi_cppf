import os
from nice_plot import *

# strain rate from 1e-2 to 1e4 per second
fdot = np.array([1e-2, 1e-1, 1.0, 10.0, 100.0, 1e3, 1e4])
# to label load file
strain_label = ['1en2', '1en1']
for i in range(2, len(fdot)):
	strain_label.append('%i' % fdot[i])

# 1st PK stress in MPa
stress = [20, 530]

def read_fp12(afile):
	# read the deformation gradient and first PK stress P of 12 from a file.
	data = np.loadtxt(afile, skiprows = 3)
	f12 = data[:, 9]
	p12 = data[:, 18]
	return f12, p12

markers = ['b-', 'g--', 'r-.', 'y:', 'm-', 'c--', 'k-.']
# fdot shown in legend
fdot_legend = [r'$10^{-2}$', r'$10^{-1}$', '$10^0$', '$10^1$', '$10^2$', '$10^3$', '$10^4$']

for j in range(len(stress)):
	plt.figure()
	# to avoid failing cases
	istart = 0
	if stress[j] == 530:
		istart = 3 # start from fdot = 10
	for i in range(istart, len(fdot)):
		print("Processing stress %i strain %s" % (stress[j], strain_label[i]))
		# back to parent dir from last submit
		os.system('cd /qfs/projects/sppsi/wkfu/grain30_grid128_increase_freq')
		# case name for label
		case = '30grain128grid_shearXY_fdot%s_p%i_power_law' % (strain_label[i], stress[j])
		fname = '30grain_128grid_30grain128grid_shearXY_fdot%s_p%i_power_law.txt' % (strain_label[i], stress[j])
		f12, p12 = read_fp12('./%s/postProc/%s' % \
							(case, fname))
		plt.plot(f12, p12 * 1e-6, markers[i], label = fdot_legend[i])
	plt.xlabel('Strain')
	plt.ylabel('Stress (MPa)')
	plt.legend(loc = 'upper left', bbox_to_anchor = (0.1, 0.9))
	plt.text(-0.22, 110, r'$\dot{F}$ (s$^{-1}$)')
	plt.savefig('ss_p%impa.png' % (stress[j]))
