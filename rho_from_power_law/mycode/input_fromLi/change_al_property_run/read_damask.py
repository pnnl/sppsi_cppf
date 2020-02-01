"""Generate gamma_dot.dat, grain_num.dat, slip_systems.dat from DAMASK outputs"""
from nice_plot import *
import os
import damask
import math
from py._builtin import enumerate
import sys
import time

# fcc 12 slip systems, copied from addSchmidfactors.py, damask source code
slip_system = np.array([
    # Slip direction     Plane normal
     [ 0, 1,-1,     1, 1, 1, ],
     [-1, 0, 1,     1, 1, 1, ],
     [ 1,-1, 0,     1, 1, 1, ],
     [ 0,-1,-1,    -1,-1, 1, ],
     [ 1, 0, 1,    -1,-1, 1, ],
     [-1, 1, 0,    -1,-1, 1, ],
     [ 0,-1, 1,     1,-1,-1, ],
     [-1, 0,-1,     1,-1,-1, ],
     [ 1, 1, 0,     1,-1,-1, ],
     [ 0, 1, 1,    -1, 1,-1, ],
     [ 1, 0,-1,    -1, 1,-1, ],
     [-1,-1, 0,    -1, 1,-1, ]])

# the cos and sin values do not change in different coordinates
# we only need the abs. values, as used in eq(4, 5) in ref1
cos_nt = np.zeros((12, 12))
cos_nb = np.zeros((12, 12))
sin_nt = np.zeros((12, 12))
sin_nb = np.zeros((12, 12))

for i in range(12):
    b_alpha = slip_system[i, 0 : 3]
    n_alpha = slip_system[i, 3 : ]
    t_alpha = np.cross(b_alpha, n_alpha)
    for j in range(12):
        b_beta = slip_system[j, 0 : 3]
        n_beta = slip_system[j, 3 : ]
        t_beta = np.cross(b_beta, n_beta)
        cos_nt[i, j] = abs(np.inner(n_alpha, t_beta) / np.linalg.norm(n_alpha) \
                    / np.linalg.norm(t_beta))
        sin_nt[i, j] = (1.0 - cos_nt[i, j] ** 2) ** 0.5
        if sin_nt[i, j] < 0.0:
            print("sin_nt[%i, %i] < 0.0" % (i, j))
        cos_nb[i, j] = abs(np.inner(n_alpha, b_beta) / np.linalg.norm(n_alpha) \
                    / np.linalg.norm(b_beta))
        sin_nb[i, j] = (1.0 - cos_nb[i, j] ** 2) ** 0.5
        if sin_nb[i, j] < 0.0:
            print("sin_nb[%i, %i] < 0.0" % (i, j))
        

class Dislocation:
    """ Calculate dislocation density from power-law results.
    ref(1): Zhao et al., IJP 80 (2016) 38-55; http://dx.doi.org/10.1016/j.ijplas.2015.12.010
    ref(2): Ma et al., Acta Materialia 54 (2006) 2169-2179; doi:10.1016/j.ijsolstr.2006.07.006
    ref(3): Ma 2004, Acta Materialia 52 (2004) 3603-3612; doi:10.1016/j.actamat.2004.04.012
    """
    
    def __init__(self, idir):
        """ idir: string 
                 folder contains damask results"""
        self.idir = idir
        self.nslip = 12 # 12 slip systems for aluminum
        
        # number of grids in three axes, int
        self.nx = 0
        self.ny = 0
        self.nz = 0
        # init nx, ny, and nz
        self.init_nxyz()
        print("nx = %i, ny = %i nz = %i" % (self.nx, self.ny, self.nz))
        self.nxyz = self.nx * self.ny * self.nz
        
        # read load file
        self.dt = 0.0
        self.init_load()
        print("dt = %.6e s" % self.dt)
        
        # init num_grain, 1D int array with size nxyz
        self.grain_map = np.zeros(self.nxyz, dtype = np.int32)
        self.init_grain_map()
        print("Grain IDs in nx = %i grids" % self.nx)
        for i in range(self.nx):
            print(self.grain_map[i])
        
        # total number of grains, normalized slip direction, plane normal, and sensor vector in the
        # global coordinate
        # self.ngrain: int
        # bsp, nsp, tsp: (self.nslip, ngrain, 3)
        self.ngrain, self.bsp, self.nsp, self.tsp = self.init_slip_systems()
        print("Number of grains = %i" % self.ngrain)
        
        # constants from Table 1, Ma 2004, Acta Materialia 52 (2004) 3603-3612
        # doi:10.1016/j.actamat.2004.04.012
        self.c1 = 0.18
        self.c2 = 5.0
        self.c3 = 5.0
        self.c4 = 8.0e6
        self.c5 = 15.0
        self.c6 = 5.0e12
        self.c7 = 1.0e-29
        self.c8 = 0.33
        self.qbulk = 2.4e-19 # in J
        self.qslip = 3.0e-19 # in J
        
        # constants to be CHECKED
        self.burger = 2.86e-10 # in m, from Phase_Nonlocal_Aluminum.config in DAMASK source code
        self.one_over_b = 1.0 / self.burger
        self.shear_modulus = 27.0e9 # in Pa, value for Al alloys, from https://www.engineeringtoolbox.com/modulus-rigidity-d_946.html
        # used to compute SSD rate, critical distance for dipole formation
        self.poisson_ratio = 0.334 # from https://www.engineeringtoolbox.com/poissons-ratio-d_1224.html
        
        self.temperature = 473.0 # in K
        self.kb = 1.38e-23 # J K-1, Boltzmann constant, from ref(3)
        # eq 24 in ref(3), or eq(9, 10) in ref(2)
        iB = 2.0 * self.kb / (self.c1 * self.c2 * self.c3 * self.shear_modulus * self.burger ** 3)
        self.b_theta =  iB * self.temperature
        
        # containers for time-dependent variables
        # The xyz grids are stored as 1D array follow DAMASK format,
        # i.e., change fastest to slowest: x, y, z
        # shear rate as a function of slip system and position
        self.gamma_dot = np.zeros((self.nslip, self.nxyz))
        # resolved shear stress as a function of slip system and position
        self.tau = np.zeros((self.nslip, self.nxyz))
        # gradient of shear rate
        self.grad_gamma_dot = np.zeros((self.nslip, self.nxyz, 3))
        
        # dislocation densities
        self.gnds = np.zeros((self.nslip, self.nxyz))
        self.gndet = np.zeros((self.nslip, self.nxyz))
        self.ssd = np.zeros((self.nslip, self.nxyz))
        self.parallel = np.zeros((self.nslip, self.nxyz))
        self.forest = np.zeros((self.nslip, self.nxyz))
        
        # dislocation density rates
        self.gnds_dot = np.zeros((self.nslip, self.nxyz))
        self.gndet_dot = np.zeros((self.nslip, self.nxyz))
        self.ssd_dot = np.zeros((self.nslip, self.nxyz))
        
    def post_process(self, incs):
        """ Post process DAMASK results
            incs: int
                  Extract results at the incs time step
        """
        
        # create folder to contain results
        if "postProc_rho" not in os.listdir(self.idir):
            mkdir = "mkdir %s/postProc_rho" % self.idir
            os.system(mkdir)
        
        # output shear rate, resolved shear stress, deformation gradient, 1st PK stress
        post_results = "postResults --time --co shearrate_slip,resolvedstress_slip --cr f,p --separation x,y,z --split " + \
                       "--range %i %i 1 -d postProc_rho *.spectralOut" % (incs, incs)
        os.system("cd %s && %s" % (self.idir, post_results))
        
        # add results to txt table
        # result txt file
        txt = "%s/postProc_rho/*inc%i.txt" % (self.idir, incs)
        
        # add gradient of shear rates
        for i in range(self.nslip):
            gradient = "addGradient --label %i_shearrate_slip %s" % ((i + 1), txt)
            os.system(gradient)
            
    def init_load(self):
        """ Read load file to get dt"""
        # find the load file
        files = os.listdir(self.idir)
        for i in files:
            if ".load" in i:
                load_file = i
                break
        load_file = self.idir + "/" + load_file
        
        # at this point, assume just one line in load file
        with open(load_file, "r") as f:
            for i, line in enumerate(f):
                if i == 0:
                    split_line = line.split()
                    for j in range(len(split_line)):
                        if split_line[j] == "time":
                            time = float(split_line[j + 1])
                        if split_line[j] == "incs":
                            incs = float(split_line[j + 1])
                            self.dt = time / incs
                            return
            
    def init_nxyz(self):
        """ Read number of grids in three axes from the geom file.
            Output:
            nx, ny, nz: int
                        number of grids in three axes
        """
        
        # find the geom file
        files = os.listdir(self.idir)
        for i in files:
            if "input_clean.geom" in i:
                geom_file = i
                break
        
        with open(self.idir + "/" + geom_file, 'r') as f:
            for i, line in enumerate(f):
                split_line = line.split()
                if split_line[0] == "grid":
                    # read number of grids in thress axes
                    # format of this line is: grid    a NX    b NY    c NZ
                    self.nx = int(split_line[2])
                    self.ny = int(split_line[4])
                    self.nz = int(split_line[6])
                    break

    def init_grain_map(self):
        """ Read from geom input.
        In damask, the microstructure indices ordered with x as fastest and z as slowest varying coordinate
        https://damask.mpie.de/Documentation/GeometryFormat
        Output:
              self.num_grain: 1D int array with size nx * ny * nz 
                              grain ID in each grid 
        """
        # find the geom file
        files = os.listdir(self.idir)
        for i in files:
            if "input_clean.geom" in i:
                geom_file = i
                break
        
        with open(self.idir + "/" + geom_file, 'r') as f:
            for i, line in enumerate(f):
                split_line = line.split()
                if i == 0:
                    # read how many lines before grain ID map
                    skip_line = int(split_line[0])
        
        # read the 2D grainID map and flatten into 1D
        self.grain_map = np.loadtxt(self.idir + "/" + geom_file, dtype = np.int32, skiprows = skip_line + 1).flatten()

    def read_column(self, txt, label):
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
        print("Located IDs " + label)
        start_time = time.time()
        data = np.loadtxt(txt, skiprows = label_line + 1, usecols = colID)
        print("time used to loadtxt = ", time.time() - start_time)
        print("end loadtxt")
        return data
    
    def read_txt(self, txt):
        """ Read the txt file generated by postResults to get shear rate, resolved shear stress, gradient shear rate."""
        with open(txt, 'r') as f:
            for i, line in enumerate(f):
                split_line = line.split()
                if split_line[0] == "inc":  # the label line
                    label_line = i
                    # the label row
                    for j in range(len(split_line)):
                        if split_line[j] == "1_shearrate_slip":
                            colID = j
                            break
                    break
        # the column arrangement is 12 shear rate, 12 resolved shear stress, and 36 gradient of shear rate
        data = np.loadtxt(txt, skiprows = label_line + 1, usecols = range(colID, colID + 12 + 12 + 36))
        for i in range(self.nslip):
            self.gamma_dot[i] = data[:, i]
            self.tau[i] = data[:, i + 12]
            for j in range(3):
                self.grad_gamma_dot[i, :, j] = data[:, 24 + i * 3 + j]
        
    def init_slip_systems(self):
        """ Read from material.config 
            Read Euler angles of grain orientation defined in material.config and transpose the vectors
            (slip direction d, plane normal n and sense vector t) in such grain-orientation coordinate to 
            reference coordinate. 
        """
        # transpose of orientation matrix, i.e., Q^T, where Q is defined in
        # Li el.al; International Journal of Plasticity 87 (2016) 154-180; Eq. (6)
        omt = []
        with open(self.idir + "/material.config", "r") as f:
            for i, line in enumerate(f):
                split_line = line.split()
                if len(split_line) > 0 and split_line[0] == "(gauss)":
                    # a sample line is:
                    # (gauss)    phi1 323.138    Phi 59.7321    phi2 342.952    scatter 0.0    fraction 1.0
                    euler_angle = [float(split_line[2]), float(split_line[4]), float(split_line[6])]
                    o = damask.Rotation.fromEulers(euler_angle, "-d")
                    omt.append(o.asMatrix())
        # number of grains
        ngrain = len(omt)
    
        # slip direction, plane normal as a function of 12 slip system and grainID. vector length 3.
        nsp = np.zeros((self.nslip, ngrain, 3))
        bsp = np.zeros((self.nslip, ngrain, 3))
        tsp = np.zeros((self.nslip, ngrain, 3))
         
        # transfer slip system into global coordinate, where vectors in local coordinate are defined in slip_system
        # v_global = O^{-1} * v_local = O^T * v_local
        for i in range(self.nslip):
            for j in range(ngrain):
                # slip direction
                b = np.dot(omt[j], slip_system[i, 0 : 3])
                b = b / np.linalg.norm(b)
                bsp[i, j] = b
                
                # plane normal
                n = np.dot(omt[j], slip_system[i, 3 : ])
                n = n / np.linalg.norm(n)
                nsp[i, j] = n
                
                # sense vector
                # according to line 2169, lattice_slip_transverse function in lattice.f90 in damask
                # source code, t = b x n, i.e., cross product.
                tsp[i, j] = np.cross(b, n)
                
        return ngrain, bsp, nsp, tsp
    
    def update_shear_rate_stress(self, txt):
        """ Read shear rate and resolved shear stress in a txt file, i.e., at a time stamp
            txt: string
                DAMASK output txt table
        """
        # loop over 12 slip systems
        for alpha in range(self.nslip):
            # shear rate and resolved stress of slip system alpha
            self.gamma_dot[alpha] = self.read_column(txt, "%i_shearrate_slip" % (alpha + 1))
            self.tau[alpha] = self.read_column(txt,"%i_resolvedstress_slip" % (alpha + 1))
    
    def update_gradient_gamma_dot(self, txt):
        """ Read the gradient of shear rate"""        
        for i in range(self.nslip):
            for j in range(3):
                label = "%i_gradFFT(%i_shearrate_slip)" % ((j + 1), (i + 1))
                self.grad_gamma_dot[i, :, j] = self.read_column(txt, label)
                
    def update_gnds_dot(self):
        """ calculate rate of GNDs dislocation density
            follow eq(7) in ref1.      
        """
        for i in range(self.nslip):
            for j in range(self.nxyz):
                grainID = self.grain_map[j]
                t = self.tsp[i, grainID - 1] # sense vector
                grad = self.grad_gamma_dot[i, j]
                value = -np.dot(grad, t) * self.one_over_b
                self.gnds_dot[i, j] = value
                
    def update_gndet_dot(self):
        """ calculate rate of GNDet dislocation density
            follow eq(7) in ref1.      
        """
        for i in range(self.nslip):
            for j in range(self.nxyz):
                grainID = self.grain_map[j]
                b = self.bsp[i, grainID - 1] # slip direction
                grad = self.grad_gamma_dot[i, j]
                value = np.dot(grad, b) * self.one_over_b
                self.gndet_dot[i, j] = value
                
    def calculate_mobile(self, parallel, forest):
        """ Calculate mobile dislocation density following eq(9) in ref2"""
        return self.b_theta * (parallel * forest) ** 0.5
    
    def calculate_tau_pass(self, rho_parallel):
        """ tau_pass = c1 * shear_modulus * burger * parallel ** 0.5
            from area between eq(3) and (4) in ref1 and
                 eq(13) in Ma 2004 by neglacting the small rho_mobile.
        """
        return self.c1 * self.shear_modulus * self.burger * rho_parallel ** 0.5
    
    def update_ssd_dot(self):
        """ calculate rate of SSN dislocation density
            follow eq(7) in ref1 and eq(24) in ref3.     
        """
        term1 = self.c4 * self.forest ** 0.5 * self.gamma_dot
        
        # term 2
        tau_pass = self.calculate_tau_pass(self.parallel)
        d_dipole = np.zeros((self.nslip, self.nxyz))
        for i in range(self.nslip):
            for j in range(self.nxyz):
                if abs(self.tau[i, j]) >= tau_pass[i, j]:
                    d_dipole[i, j] = 3.0 ** 0.5 * self.shear_modulus * self.burger \
                    / 16.0 / np.pi / (1.0 - self.poisson_ratio) / self.tau[i, j]
        
        mobile = self.calculate_mobile(self.parallel, self.forest)
        term2 = self.c6 * d_dipole * mobile * self.gamma_dot
        
        term3 = -self.c5 * self.ssd * self.gamma_dot
        
        # const of term 4
        consts = -self.c7 * \
        math.exp(-self.qbulk / self.kb / self.temperature) / \
        self.kb / self.temperature
        # von Mises equi. stress = absolute value of tau, eq (19), ref2
        # von Mises equi. shear rate using the abs. value at this point
        # TO DO: check von Mises equi. definition of shear rate
        term4 = consts * np.absolute(self.tau) * self.ssd ** 2 * \
        np.absolute(self.gamma_dot) ** self.c8
        
        self.ssd_dot = term1 + term2 + term3 + term4
        
    def update_forest(self):
        """ Update forest dislocation density following eq(4) in ref1
        """
        # set to zero
        self.forest.fill(0.0)
        for i in range(self.nslip):
            for j in range(self.nslip):
                # note the absolute cosines are stored
                self.forest[i] += self.ssd[j] * cos_nt[i, j] + \
                                  np.absolute(self.gnds[j] * cos_nb[i, j]) + \
                                  np.absolute(self.gndet[j] * cos_nt[i, j])
    
    def update_parallel(self):
        """ Update parallel dislocation density following eq(5) in ref1
        """
        # set to zero
        self.parallel.fill(0.0)
        for i in range(self.nslip):
            for j in range(self.nslip):
                self.parallel[i] += self.ssd[j] * sin_nt[i, j] + \
                                  np.absolute(self.gnds[j] * sin_nb[i, j]) + \
                                  np.absolute(self.gndet[j] * sin_nt[i, j])
    
    def time_integrator(self, v0, rate, dt):
        """ Compute v1 as
              v1 - v0 = dt * rate(at t1), eq(50) in damask review paper
            v0 and rate must be in the same dimension.
            rate should be or close to t1, the updated value 
        """
        v1 = v0 + rate * dt
        return v1
    
    def sum_over_12_slip_system(self, rho):
        """ Sum over slip systems"""
        s_rho = np.zeros(self.nxyz)
        for i in range(self.nslip):
            s_rho += rho[i]
        return s_rho   
    
    def solver(self):
        """ Calculate dislocation density evolution.
        """
        # form the txt files in the sequence of increment
        # the three processed steps for testing
        files = os.listdir(self.idir + "/postProc_rho/")
        
        # lists to store sum of gnd and ssd
        gnd_list = []
        ssd_list= []
        
        for i in range(1, 101):
            print("Process incs %i" % i)
            time_inc = time.time()
            # locate txt file
            for j in files:
                if "law_inc%i.txt" % i in j:
                    txt = self.idir + "/postProc_rho/" + j
                    break
            
            # update shear rate, resolved shear stress and gradient of shear rate
            if i == 1:
                print("start read txt")
                start_time = time.time()
            self.read_txt(txt)
            if i == 1:
                print("end read txt")
                print("used time = ", time.time() - start_time)

            # update GNDs and GNDet rate
            self.update_gnds_dot()
            self.update_gndet_dot()
            
            # calculate GNDs and GNDet
            self.gnds = self.time_integrator(self.gnds, self.gnds_dot, self.dt)
            self.gndet = self.time_integrator(self.gndet, self.gndet_dot, self.dt)
            
            gnd = (self.gnds ** 2 + self.gndet ** 2) ** .5
            
            # sum over 12 slip systems
            gnd_total = self.sum_over_12_slip_system(gnd)
            np.savetxt(self.idir + "/postProc_rho/gnd_inc%i.txt" % i, gnd_total)
            sum_gnd = np.sum(gnd_total)
            gnd_list.append(sum_gnd)
            print("Sum of GND is %.6e" % sum_gnd)
            np.savetxt(self.idir + "/postProc_rho/gnd_sum.txt", gnd_list)
            
            # use the updated GNDs and GNDet to update parallel and forest
            # i.e., parallel and forest at half step
            self.update_forest()
            self.update_parallel()
            
            # TO DO: since ssd_dot depends on parallel and forest, an iterative scheme might be
            # needed  
            # use the opdated, half-step parallel and forest to compute ssd_dot
            self.update_ssd_dot()
            
            # update ssd
            self.ssd = np.absolute(self.time_integrator(self.ssd, self.ssd_dot, self.dt))
            ssd_total = self.sum_over_12_slip_system(self.ssd)
            np.savetxt(self.idir + "/postProc_rho/ssd_inc%i.txt" % i, ssd_total)
            ssd_list.append(np.sum(ssd_total))
            np.savetxt(self.idir + "/postProc_rho/ssd_sum.txt", ssd_list)

            # update parallel and forest
            self.update_forest()
            np.savetxt(self.idir + "/postProc_rho/forest_inc%i.txt" % i, self.sum_over_12_slip_system(self.forest))
            self.update_parallel()
            np.savetxt(self.idir + "/postProc_rho/parallel_inc%i.txt" % i, self.sum_over_12_slip_system(self.parallel))
            
            time_inc = time.time() - time_inc
            print("Inc %i used %2f s" % (i, time_inc))

# dirname = "./power/postProc_rho"
# for i in range(3):
#     forest = np.loadtxt(dirname + "/forest_inc%i.txt" % i)
#     parallel = np.loadtxt(dirname + "/parallel_inc%i.txt" % i)
#     
#     # sum over 12 slip systems
#     sforest = np.zeros(len(forest[0]))
#     sparallel = np.zeros(len(parallel[0]))
#     for j in range(len(forest)):
#         sforest += forest[j]
#         sparallel += parallel[j]
#         
#     plt.figure()
#     plt.imshow(np.reshape(sforest[0 : 16 * 16], (16, 16)), aspect = 'auto')
#     plt.colorbar()
#     plt.xlabel("x")
#     plt.ylabel("y")
#     plt.savefig("forest_inc%i.png" % i)
#     
#     plt.figure()
#     plt.imshow(np.reshape(sparallel[0 : 16 * 16], (16, 16)), aspect = 'auto')
#     plt.colorbar()
#     plt.xlabel("x")
#     plt.ylabel("y")
#     plt.savefig("parallel_inc%i.png" % i)



