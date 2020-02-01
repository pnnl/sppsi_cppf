"""Test the euler angle definition in DAMASK
literature: Li el.al; International Journal of Plasticity 87 (2016) 154-180; Eq. (6) 
"""
import damask
import numpy as np

# phi1, Phi, phi2 in degree and change to radians
angle_deg = np.array([100.0, 20.0, 140.0])
angle_rad = angle_deg * np.pi / 180.0

c = np.cos(angle_rad)
s = np.sin(angle_rad)

om = np.zeros((3, 3))

# codes copy from eu2om(eu) defined in rotations.f90, damask
om[1-1, 1-1] = c[1-1]*c[3-1]-s[1-1]*s[3-1]*c[2-1]
om[1-1,2-1] =  s[1-1]*c[3-1]+c[1-1]*s[3-1]*c[2-1]
om[1-1,3-1] =  s[3-1]*s[2-1]
om[2-1,1-1] = -c[1-1]*s[3-1]-s[1-1]*c[3-1]*c[2-1]
om[2-1,2-1] = -s[1-1]*s[3-1]+c[1-1]*c[3-1]*c[2-1]
om[2-1,3-1] =  c[3-1]*s[2-1]
om[3-1,1-1] =  s[1-1]*s[2-1]
om[3-1,2-1] = -c[1-1]*s[2-1]
om[3-1,3-1] =  c[2-1]

print("om using equations")
print(om)

o = damask.Rotation.fromEulers(angle_deg, "-d")

print("om using damask module")
print(o.asMatrix())

# conclusion: the output of the above command is the transpose of Q defined in literature (line 2)
# can be directly used to calculate vector defined in local coordinate to global coordinate