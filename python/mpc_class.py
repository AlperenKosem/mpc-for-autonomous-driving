import numpy as np
import scipy
from scipy import signal
import matplotlib.pyplot as plt
from qpsolvers import solve_qp


class mpc_controller :
    def __init__(self, Nc, Np):

        self.trajectory_points = None 
        self.transform = None
        self.vehicle_state = None
        self.rw = 0.1 # input weight
        self.amplitude_constraint = None
        self.rate_constraint = None

        # vehicle parameters

        self.Cf = 155494.663 
        self.Cr = 155494.663
        self.Iz = 1436.24
        self.Lf = 1.165
        self.Lr = 1.165
        self.m = 1140
        self.Vx = 30
        self.dt = 0.01 # sample time

        self.Np = Np   # Predict Horizon
        self.Nc = Nc   # Control Horizon
        
        # DYNAMIC MODEL IN TERMS OF ERROR WITH RESPECT TO ROAD
       
        #  states : distance of the line center, distance rate of line center,
        #  orientation error, orientation error rate
        #  inputs : desired yaw rate(Vx/R), steering angle
        self.A = np.zeros((4,4))
        self.A[0,1] = 1
        self.A[1,1] = -(self.Cf + self.Cr) / (self.m * self.Vx)
        self.A[1,2] = ( self.Cf +  self.Cr) / self.m
        self.A[1,3] = (- self.Cf * self.Lf +  self.Cr * self.Lr) / (self.m * self.Vx)
        self.A[2,3] = 1
        self.A[3,1] = -( self.Cf * self.Lf -  self.Cr * self.Lr) / (self.Iz * self.Vx)
        self.A[3,2] = ( self.Cf * self.Lf - self.Cr * self.Lr) / self.Iz
        self.A[3,3] = -( self.Cf * self.Lf**2 +  self.Cr * self.Lr**2) / (self.Iz * self.Vx)

        
        self.B = np.zeros((4,2))
        self.B[1,0] = ( self.Cf/self.m) 
        self.B[1,1] = (-( self.Cf * self.Lf -  self.Cr * self.Lr)/(self.m * self.Vx)) - self.Vx
        self.B[3,0] = ( self.Cf * self.Lf) / self.Iz
        self.B[3,1] = -( self.Cf * self.Lf**2 +  self.Cr * self.Lr**2)/ (self.Iz * self.Vx)

        self.C = np.zeros((2,4))
        self.C[0,0] = 1
        self.C[1,2] = 1

        self.D = np.zeros((2,2)) 

        # A B C D Checked.

        self.nstates = self.A.shape[0]  # number of states
        self.ninputs = None             # number of inputs
        self.noutputs = self.C.shape[0] # number of outputs

        ## Full State Feedback

        self.P = np.array([-5-3j, -5+3j , -7 , -10 ])
        self.B1 = self.B[:,0]
        self.B1 = np.reshape(self.B1,(self.nstates,1))
        fsf1 = signal.place_poles(self.A, self.B1, self.P,)
        K = fsf1.gain_matrix
        self.A = (self.A - self.B1 * K) # state feedback controller A

        self.discrete_system = scipy.signal.cont2discrete((self.A,self.B,self.C,self.D), self.dt, method='zoh')
        
        self.Ad = self.discrete_system[0]
        self.Bd = self.discrete_system[1]
        
                
        self.Bd1 = self.Bd[:,0]
        self.Bd2 = self.Bd[:,1]
        self.Bd1 = np.reshape(self.Bd1,(self.nstates,1))
        self.Bd2 = np.reshape(self.Bd2,(self.nstates,1))
        self.Cd = self.discrete_system[2]
        self.Dd = self.discrete_system[3]

        self.ninputs = len(self.Bd1[0]) 
    
    
    def generate_augmented_model(self):

        ##### A Augmented
        self.A_aug = np.zeros((self.nstates + self.noutputs, self.nstates + self.noutputs))

        for i in range(0,self.nstates)  :
            for j in range(0,self.nstates)  :
                self.A_aug[i,j] = self.Ad[i,j]
        
        CdxAd = np.matmul(self.Cd, self.Ad)


        for i in range(0,self.noutputs) :
            for j in range(0,self.nstates) :
                
                self.A_aug[self.nstates + i,j] = CdxAd[i,j]
    
        for i in range(self.nstates,self.nstates + self.noutputs) :
            for j in range(self.nstates, self.nstates + self.noutputs) :
                if i == j  :
                    self.A_aug[i,j] = 1
                
        
        #### B Augmented

        self.B_aug = np.zeros((self.nstates + self.noutputs, self.ninputs))
        CdxBd = np.matmul(self.Cd, self.Bd1)
        for i in range(0,self.nstates) :
            for j in range(0,self.ninputs) :
                    self.B_aug[i,j] = self.Bd1[i,j]
                  

        for i in range(0, self.noutputs) :
            for j in range(0,  self.ninputs) :
                self.B_aug[self.nstates + i,j] = CdxBd[i,j]
                # self.B_aug[i,j] = CdxBd[i - self.nstates,j - self.nstates]

        #### Bw Augmented

        self.Bw_aug = np.zeros((self.nstates + self.noutputs, self.ninputs))
        CdxBdw = np.matmul(self.Cd, self.Bd2)
        for i in range(0,self.nstates) :
            for j in range(0,self.ninputs) :
                    self.Bw_aug[i,j] = self.Bd2[i,j]
                  

        for i in range(0, self.noutputs) :
            for j in range(0,  self.ninputs) :
                self.Bw_aug[self.nstates + i,j] = CdxBdw[i,j]


        #### C Augmented
        
        self.C_aug = np.zeros((self.noutputs,self.nstates + self.noutputs))
        # self.C_aug[0, self.nstates ] = 1

        for i in range(0, self.noutputs) :
            for j in range(self.nstates,self.nstates + self.noutputs) :
                if j - i == self.nstates  :
                    self.C_aug[i,j] = 1

        # D Augmented
        self.D_aug = [0] #unused
    
        # returns self.A_aug, self.B_aug, self.C_aug, self.D_aug

    def calculate_mpc_gains(self):

        C_augxA_aug = np.matmul(self.C_aug, self.A_aug) # for calculation F matrix
        rows = len(C_augxA_aug)
        cols = len(C_augxA_aug[0])

        self.F = np.zeros((rows * self.Np, cols)) 
        pow = 1
        m = 0
        j = 0
        
        for i in range(0,self.noutputs * self.Np, 2 ) :

            dummy = np.matmul(self.C_aug, np.linalg.matrix_power(self.A_aug,pow))
            
            while j  <  cols :
                for k in range(0,self.noutputs) :
                    self.F[i + k ,j ] = dummy[k,j]
                j += 1
            j = 0
            pow +=1

        self.Phi = np.zeros((self.Np  * self.noutputs, self.Nc * self.ninputs))    # phi created for input
        self.Omega = np.zeros((self.Np  * self.noutputs, self.Nc * self.ninputs))  # self omega created for road curvature
        i = 0
        k = 0
        m = 0
        j = 0
        for j in range(0,self.ninputs * self.Nc):
            while( i < self.Np) :
                
                pow_a = np.linalg.matrix_power(self.A_aug,i - k)
                dummy1 = np.matmul(self.C_aug, pow_a)
                dum_b = np.matmul(dummy1,self.B_aug)
                dum_bw = np.matmul(dummy1, self.Bw_aug)

                self.Phi[m, j] = dum_b[0]
                self.Phi[m + 1, j] = dum_b[1]

                self.Omega[m, j] = dum_bw[0]
                self.Omega[m + 1, j] = dum_bw[1]
                i += 1
                m += 2

            i = j + 1
            m = j + i + 1
            k = k + 1 

        self.PhiT_phi = np.matmul(np.transpose(self.Phi), self.Phi)
        self.PhiT_F = np.matmul(np.transpose(self.Phi), self.F)
        self.PhiT_Omega = np.matmul(np.transpose(self.Phi), self.Omega)
        self.PhiT_Rs = self.PhiT_F[:, -self.noutputs:] 
    
    def quad_solver(self, ref_act, Xf, deltaW) : 

        ## set constraints
        self.M = np.zeros((self.Nc * 4, self.Nc), dtype=np.double)
        self.gamma = np.zeros((self.Nc * 4,1), dtype=np.double)
        for i in range(0,self.Nc) :
            for j in range(0, self.Np) :
                if i == j :
                    self.M[i,j] = 1
                    self.gamma[i,] = self.rate_constraint
                    self.M[i + self.Nc,j ] = -1
                    self.gamma[i + self.Nc, ] = self.rate_constraint
                if i >= j :
                    self.M[i + 2 * self.Nc, j] = 1
                    self.M[i + 3 * self.Nc, j] = -1
                    self.gamma[ i + 2 * self.Nc, ] = self.amplitude_constraint - self.u0
                    self.gamma[ i + 3 * self.Nc, ] = self.amplitude_constraint + self.u0
        
        E = 2 * (self.PhiT_phi + self.rw * np.eye(self.Nc))
        Xf = np.reshape(Xf, (6,))
        deltaW = np.reshape(deltaW, (self.Nc,))
        Fconstraint = -2 * (np.matmul(self.PhiT_Rs , ref_act) - np.matmul(self.PhiT_F , Xf) - np.matmul(self.PhiT_Omega , deltaW ))
        
        self.gamma = np.reshape(self.gamma, (self.Nc * 4,))
        Fconstraint = np.reshape(Fconstraint, (self.Nc,))
        x = solve_qp(E, Fconstraint, self.M, self.gamma)
        return x[0]
        



    def mpc_simulation(self,N_sim, ref_signal, curvature, u0, amplitude_constraint, rate_constraint, xm , xf) :
        self.u0 =  u0
        # self.u_act = 
        self.amplitude_constraint = amplitude_constraint
        self.rate_constraint = rate_constraint



        u2 = np.ones((1, N_sim + Nc))  * (curvature * self.Vx)
        deltaW = np.zeros((Nc,1))
        self.u = np.zeros_like(u2)
        self.deltaUs = np.zeros_like(u2)
        self.y = np.zeros((2, N_sim))
        self.y[0,0] = xm[0,0]
        self.y[1,0] = xm[2,0]
        for i in range(1,N_sim - 1) :
            
            ref_act = ref_signal[:, i]
            self.deltaU = self.quad_solver(ref_act, xf, deltaW)
                       
            self.deltaUs[0,i] = self.deltaU
            self.u[0,i] = self.u[0, i-1] + self.deltaU
            
            xm_old = xm
            xm = np.matmul(self.Ad, xm)  + self.Bd1 * self.u[0,i] + self.Bd2 * u2[0,i]
            ans1 = np.matmul(self.Cd , xm)

            self.y[0,i] = ans1[0]
            self.y[1,i] = ans1[1]
            a = xm - xm_old
            xf[0] = a[0]
            xf[1] = a[1]
            xf[2] = a[2]
            xf[3] = a[3]
            xf[4] = self.y[0,i]
            xf[5] = self.y[1,i]
            self.u0 = self.u[0,i]
            
        
if __name__ == '__main__':

    # Parameters for sim
    Np = 12
    Nc = 3
    
    N_sim = 300
    curvature = 1/ 100
    ref_signal = np.zeros((2, N_sim)) 
    xm = np.zeros((4, 1)) 
    xf = np.zeros((6,1))
    xm[0,0] = 1.5
    xf[4,0] = 1.5
    xm[2,0] = np.radians(-30)
    xf[5,0] = np.radians(-30)
    
    u0 = [0] 
    amplitude_constraint = np.radians(35)
    rate_constraint = np.radians(15) 

    controller = mpc_controller(Nc,Np)
    controller.generate_augmented_model()
    controller.calculate_mpc_gains()

    controller.mpc_simulation(N_sim, ref_signal, curvature,u0 ,amplitude_constraint, rate_constraint, xm, xf )
    
    ## plotting 
    plt.figure(1)
    len_u = len(np.reshape(controller.u,(N_sim + Nc,1)))
    plt.plot(np.reshape(controller.u,(N_sim + Nc,1)))
    plt.plot(amplitude_constraint * np.ones((len_u,1)), 'r--')
    plt.plot( - amplitude_constraint * np.ones((len_u,1)), 'r--')
    plt.ylabel('steering angle( radian )')
    plt.xlabel(' sample ')
    plt.grid(True)

    plt.figure(2)
    len_u = len(np.reshape(controller.deltaUs,(N_sim + Nc,1)))
    plt.plot(np.reshape(controller.deltaUs,(N_sim + Nc,1)))
    plt.plot(rate_constraint * np.ones((len_u,1)), 'r--')
    plt.plot( - rate_constraint * np.ones((len_u,1)), 'r--')
    plt.ylabel('steering angle( radian )')
    plt.xlabel(' sample ')
    plt.grid(True)
    

    plt.figure(3)
    plt.plot(np.reshape(controller.y[0,:],(N_sim,1)), label = 'system response')
    plt.plot(ref_signal[0,:], label = 'reference signal')
    plt.legend()
    plt.ylabel('distance of lane center')
    plt.xlabel(' sampling ')
    plt.grid(True)
 

    plt.figure(4)   
    plt.plot(np.reshape(controller.y[1,:],(N_sim,1)), label = 'system response')
    plt.plot(ref_signal[1,:],  label = 'reference signal' )
    plt.legend()
    plt.ylabel(' yaw angle respect to lane center ( radian )')
    plt.xlabel(' sampling ')
    plt.grid(True)

    plt.show()


