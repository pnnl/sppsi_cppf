import os
""" Generator vtr file to show dislocation density distribution in Paraview"""

class Paraview_generator:
    def __init__(self, idir):
        self.idir = idir
        
    def generator_vtr_container(self, inc):
        folder = self.dir + "/postProc_rho"
        files = os.listdir(folder)
        for file in files:
            if "law_inc%i.txt" % inc in file:
                txt = file
                break
        cmd = "vtk_rectilinearGrid %s" % txt
        os.system(cmd)