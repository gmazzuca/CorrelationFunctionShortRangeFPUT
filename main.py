''' Performe the evolution of the FPUT chain with potential of the form:
H = \sum_{j=0}^{N-1} (p_j^2/2 + \sum_{s=1}^m k_s((q_{j+s} - q_j)^2/2 + chi*(q_{j+s} - q_j)^3/3 + gamma*(q_{j+s} - q_j)^4/4) ))
We will refer to k_s as the springs intensities.
The evolution is made with a yoshida algorithm of order 4 realize through a Leap Frog.
Note that the code is tought to be used on a cluster running python3.
To work the program needs the following arguments:
number = int, denote the number of file
beta = double, inverse of temperature
tfinale = double, final time of the evolution
particles= int, number of particles
trials = int, number of trials
chi = double parameter in front of the cubic part of the potential
gamma = double, parameter in front of the quortic part of the potential
integration_step = double (small), integration step
log_step = double, logarithmic spacing, the code saves the values of the correlation functions  in logarithmic scale
'''

import numpy as np
import sys



################### MISCELLANEA ####################################################
''' Not specific functions for the evolution, but useful in order to sample initial data '''
def DHT(x):
    ''' Discrete Hartley Transform'''
    fx = np.fft.fft(x)
    return np.real(fx) - np.imag(fx)


def DHTn(x):
    ''' Discrete Hartley Transform NORMALIZED'''
    fx = np.fft.fft(x)
    return (np.real(fx) - np.imag(fx))/np.sqrt(len(x))


def space_log(t,logstep):
    ''' logarithmic scale of time '''
    
    xfine = np.arange(0,int(min([t,20])))
    while (xfine[-1] < t):
        tmp = xfine[-1]*logstep
        if tmp < t:
            xfine = np.append(xfine, tmp)
        else:
            xfine = np.append(xfine,t)

    return xfine


def periodic_difference(x,d):
    ''' give the periodic distance at d level '''
    ytmp = np.append(x[-d:],np.append(x,x[:d]))
    y = ytmp[d:] - ytmp[:-d]
    return y[d:]

def shift(A,d):
    ''' shifting elements of A's row by a their number '''
    B = np.zeros(A.shape)
    for l in range(d):
        B[l,:] = np.append(A[l,-l-1:], A[l,:-l-1])

    return B
############# Check functions (energy) #############


def periodic_energy(p,q,d,coef):

    '''energy periodic toda in (p,q) variables'''
    A = np.zeros((d,len(q)))
    Acoef  = np.zeros((d,len(q)))
    for l in range(d):
        A[l,:] =  periodic_difference(q,l+1)
        Acoef[l,:] = coef[l]*A[l,:]
    qpart = np.sum( A*Acoef*(0.5 + A*(chi/3 + A*gamma*0.25)))
    ppart = 0.5*np.sum(p*p)
    
    return qpart + ppart


######### Evolution with a LP algorithm ###############
def vecp(q,d,coef):
    
    ''' force for the evolution, must give dH/dq '''
    A = np.zeros((d,len(q)))
    Acoef  = np.zeros((d,len(q)))
    for l in range(d):
        A[l,:] =  periodic_difference(q,l+1)
        Acoef[l,:] = coef[l]*A[l,:]
    forcematrix =  Acoef*(1 + A*(chi + A*gamma))
    force = np.sum(- forcematrix + shift(forcematrix,d) ,0)
    return force

def leap_frog(p,q,dt,d,coef):
    ''' Leap frog generic '''
    
    q_tmp = q + 0.5*dt*p
    p_new = p -  dt*vecp(q_tmp,d,coef)
    q_new = q_tmp + 0.5*dt*p_new

    return(p_new,q_new)

def yo4(p,q,dt,d,coef):
    ''' time step integration via yoshida4'''

    x1= 1.351207191959657
    x0 = -1.702414383919315
    
    (p1, q1) = leap_frog(p,q, x1*dt,d,coef)
    (p2, q2) = leap_frog(p1,q1, x0*dt,d,coef)
    (p3, q3) = leap_frog(p2,q2, x1*dt,d,coef)

    return (p3,q3)

def evolution_LP(p,q,tau,dt,d,coef):
    ''' periodic evolution till time tau of the data (p,q) with a timestep dt performed with a Leap Frog algorithm of order 2'''
    time = 0
    while(time < tau):
        (p,q) = leap_frog(p,q,dt,d,coef)
        time = time + dt
        
    return(p,q)

def evolution_yo4(p,q,tau,dt,d,coef):
    ''' periodic evolution till time tau of the data (p,q) with a timestep  dt performed with a Yoshida algorithm of order 4'''
    time = 0
    while(time < tau):
        (p,q) = yo4(p,q,dt,d,coef)
        time = time + dt
        
    return(p,q)

def complete_sol(p,q,time,dt,d,coef ):
    ''' complete solution for the periodic case '''
    tsteps = len(time)
    particles = len(p)

    solp = np.zeros((tsteps, particles))
    solq = np.zeros((tsteps, particles))
    evostep = np.ediff1d(time)

    solp[0,:] = p
    solq[0,:] = q
    
    for k in range(tsteps - 1 ):
        (solp[k+1],solq[k+1]) = evolution_yo4(solp[k],solq[k],evostep[k], dt,d,coef)

    return (solp,solq)


def circmatrix(a):
    ''' generate a circulant matrix starting from vector '''
    n = len(a)

    M = np.zeros((n,n))
    for k in range(n):
        M[k,:] = np.append(a[k:], a[:k])

    return M

def circ_root(eig):

    ''' square root of circulant matrix '''
    for k in range(len(eig)):
        if eig[k] < 0 :
            print('impossibile fare la radice\n')
            exit()


    sqrteig = np.sqrt(eig)
    vector = DHT(sqrteig)/len(eig)

    return circmatrix(vector)



######### Starting Point - Stats #####################


def eigenvalues_q(coef,n):
    ''' eigenvalues of the interacting matrix of q '''

    d = len(coef)
    compl_coef = np.append(np.append(np.append(2*sum(coef), -coef) , np.zeros(n - 2*d-1)) , - coef[::-1])
    eig = DHT(compl_coef)
    eig[0] = 0
    return eig


def initial_condition(n,beta,eig):

    sigma = 1/np.sqrt(beta)

    # INITIAL CONDITION ON P
    tmp_p = np.random.normal(loc=0.0, scale=sigma, size= n-1) # normal on the independent variable
    tmp_p = np.append(0,tmp_p)
    p = DHTn(tmp_p)

    # INITIAL CONDITION ON Q

    tmp_q = np.random.normal(loc=0.0, scale=sigma/np.sqrt(eig[1:]), size= n-1) # normal on the independent variable
    tmp_q = np.append(0,tmp_q)
    q = DHTn(tmp_q)

    return (p,q)

    
    
################# MAIN ##########################à


if (len(sys.argv) < 9):
    print('error: give number beta tfinale particles trials chi gamma integration_step log_step\n')
    exit()


global chi
global gamma

np.random.seed()
# PARAMETERS FOR EVOLUTION FROM COMMAND LINE
number = int(sys.argv[1])
beta = int(sys.argv[2])
tfinal = float(sys.argv[3])
n = int(sys.argv[4])
trials = int(sys.argv[5])
chi = float(sys.argv[6])
gamma = float(sys.argv[7])
dt = float(sys.argv[8])
logstep = float(sys.argv[9])

# springs intensities vector
coef = np.array([1,0.125, 7/72])
d = len(coef)

# setting time e strobes
timer = space_log(tfinal,logstep)

numstrobe = len(timer)
    
# Setting Storage Vector

base_index = int(np.ceil(n*0.5)) #central particles
#correlation
#S22 refers to che correlation Cor(p_j(t), p_0(0))
#S11 refers to che correlation Cor(r_j(t), r_0(0))
#S12 refers to che correlation Cor(r_j(t), p_0(0))
#S21 refers to che correlation Cor(p_j(t), r_0(0))
S22 = np.zeros((numstrobe,n))
S11 = np.zeros((numstrobe,n))
S12 = np.zeros((numstrobe,n))
S21 = np.zeros((numstrobe,n))
r_mean_values_evolution = np.zeros((numstrobe,n)) #contains all the mean values
p_mean_values_evolution = np.zeros((numstrobe,n))



############ Evolution ##################
eig_q_matrix = eigenvalues_q(coef,n) # eigenvalues for the matrix
M = circ_root(eig_q_matrix) 
for k in range(n):
    if eig_q_matrix[k] < 0 :
        print('error!! negative eigenvalue')
        exit()

################## Energy, uncomment if you need to check that the code works properly ##########
#e = np.zeros(numstrobe)
# evolutions 
for k in range(trials):
    
    (p,q) = initial_condition(n,beta,eig_q_matrix)
    p0 = p[base_index]
    r0 = np.dot(M[base_index,:],q)
    (psol,qsol) = complete_sol(p,q,timer,dt,d,coef)
    
    # Energy check, uncomment if you need it
    #for l in range(numstrobe):
     #   e[l] = periodic_energy(psol[l,:],qsol[l,:],d,coef)
    #plt.plot(timer, (e-e[0])/e[0],'b')
    #plt.show()
     #statistical quantities
    rsol = np.matmul(qsol,M)
    r_mean_values_evolution = r_mean_values_evolution + rsol  #to evaluate the mean value of single particles
    p_mean_values_evolution = p_mean_values_evolution + psol
    S22 = S22 +  psol * p0 # main part of correlation
    S11 = S11 + rsol * r0
    S12 = S12 + rsol*p0
    S21 = S21 + psol*r0
    


# Computing the averages:

r_mean_values_evolution = r_mean_values_evolution/trials
p_mean_values_evolution = p_mean_values_evolution/trials
S22 = S22/trials
S11 = S11/trials
S21 = S21/trials
S12 = S12/trials

S22 = S22 - p_mean_values_evolution*p_mean_values_evolution[0,base_index]
S11 = S11 - r_mean_values_evolution*r_mean_values_evolution[0,base_index]
S21 = S21 - p_mean_values_evolution*r_mean_values_evolution[0,base_index]
S12 = S12 - r_mean_values_evolution*p_mean_values_evolution[0,base_index]




# Saving Files
np.savetxt('FPU_S22_n_%d_beta_%d_chi_%f_gamma_%f_%05d.dat' %(n,beta,chi,gamma,number) ,S22)
np.savetxt('FPU_S11_n_%d_beta_%d_chi_%f_gamma_%f_%05d.dat' %(n,beta,chi,gamma,number) ,S11)
np.savetxt('FPU_S12_n_%d_beta_%d_chi_%f_gamma_%f_%05d.dat' %(n,beta,chi,gamma,number) ,S12)
np.savetxt('FPU_S21_n_%d_beta_%d_chi_%f_gamma_%f_%05d.dat' %(n,beta,chi,gamma,number) ,S21)



