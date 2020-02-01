import os
import numpy as np

def mat_generator():
    s = """#-------------------#
<homogenization>
#-------------------#
[directSX]
mech    none

#-------------------#
<crystallite>
#-------------------#
[almostAll]
(output) phase
(output) texture
(output) volume
(output) orientation    # quaternion
(output) grainrotation  # deviation from initial orientation as axis (1-3) and angle in degree (4)
(output) f              # deformation gradient tensor; synonyms: "defgrad"
(output) fe             # elastic deformation gradient tensor
(output) fp             # plastic deformation gradient tensor
(output) p              # first Piola-Kichhoff stress tensor; synonyms: "firstpiola", "1stpiola"
(output) lp             # plastic velocity gradient tensor

<phase>
{Phase_Phenopowerlaw_Aluminum.config}

<texture>
"""

    # read euler angle of 40 grains
    angle = np.loadtxt("euler_angle.dat")
    
    for i in range(len(angle)):
        # in DAMASK, phi1 and phi2 from 0 to 2pi, Phi from 0 to pi
        # in euler_angle.dat, phi1 and phi2 from -pi to pi, Phi from -pi/2 to pi/2
        # to be consistent
        # if phi < 0, phi = 360. + phi
        # if Phi < 0, Phi = -Phi + 90
        phi1 = angle[i][1]
        if phi1 < 0.0:
            phi1 = 360.0 + phi1
        phi2 = angle[i][3]
        if phi2 < 0.0:
            phi2 = 360.0 + phi2
        Phi = angle[i][2]
        if Phi < 0.0:
            Phi = -Phi + 90.0
        s += """[Grain%i]
(gauss) phi1 %.2f Phi %.2f phi2 %.2f scatter 0.0 fraction 1.0\n""" % ((i + 1), phi1, Phi, phi2)
        
    s += "\n<microstructure>\n"
        
    for i in range(len(angle)):
        s += """[Grain%i]
crystallite 1
(constituent)    phase 1    texture %i    fraction 1.0
""" % ((i + 1), (i + 1))
    
    with open("material.config", "w") as f:
        f.write(s)
        
def geom_generator():
    id = np.loadtxt("grain_num.dat")
    print("length is ", len(id))
    # length along one axis
    ngrid = 64
    s = """4 header
grid a 64 b 64 c 64
size x 1  y 1  z 1
origin x 0.0 y 0.0 z 0.0
homogenization 1
"""
    for i in range(len(id)):
        s += "%i " % id[i]
        if i % ngrid == 0 and i > 0:
            s += "\n"
    
    with open("input.geom", "w") as f:
        f.write(s)
    

mat_generator()
geom_generator()







    