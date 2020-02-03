from nice_plot import *
import os

""" Post process DAMASK results to plot stress-strain curve"""

class SS_curve:
    def __init__(self, idir):
        """ idir: string 
                  folder contains the spectralOut file"""
        self.idir = idir
        
    def add_ss(self):
        """ Process the txt files to add Cauchy stress, strain, and their von Mises."""
        # The txt file is in postProc folder in self.idir
        folder = self.idir + "/postProc"
        assert len(os.listdir(folder)) == 1
        txt = os.listdir(folder)[0]
        cmd1 = "addCauchy %s" % txt
        cmd2 = "addStrainTensors --left --logarithmic %s" % txt
        cmd3 = "addMises -s Cauchy %s" % txt
        cmd4 = 'addMises -e "ln(V)" %s' % txt
        os.system("cd %s && %s && %s && %s && %s" % (folder, cmd1, cmd2, cmd3, cmd4))

# dirs = ["30grain128grid_shearXY_fdot1000_p20_power_law", "30grain128grid_shearXY_fdot1_p20_power_law", "30grain128grid_shearXY_fdot1en1_p20_power_law", 
#         "30grain128grid_shearXY_fdot1en2_p20_power_law"]
# for idir in dirs:
#     a = SS_curve(idir)
#     a.add_ss()