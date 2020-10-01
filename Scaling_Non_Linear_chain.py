''' This file is used to make the plot of the correlations functions. It reads the file that FPUT_short_range.py produce.
We recall that we are studing the correlation function of system of the form
H = \sum_{j=0}^{N-1} (p_j^2/2 + \sum_{s=1}^m k_s((q_{j+s} - q_j)^2/2 + chi*(q_{j+s} - q_j)^3/3 + gamma*(q_{j+s} - q_j)^4/4) ))

Parameters from command line:

beta = int, inverse of the temperature
tfinal = float, final time of the evolution
n = int number of particles
chi = float, parameter in front of the cubic part of the potential 
gamma = float, parameter in front of the quartic part of the potential
logstep = float, logarithmic step, it is used to save data

'''
import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
import scipy.optimize as opt
import matplotlib.pylab as pylab
params = {'legend.fontsize': 30,
          'figure.figsize': (19, 13),
         'axes.labelsize': 45,
         'axes.titlesize':50,
         'xtick.labelsize':25,
          'ytick.labelsize':25}
pylab.rcParams.update(params)



##### Compute how the max of correlation scales logaritmic and its speed

def func(x,a,b):
    return b*x**a

def speed(x,a):
    return a*x

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

if (len(sys.argv) < 8):
    print('error: give  beta tfinale particles  chi gamma log_step example_number')
    exit()

# PARAMETERS FOR EVOLUTION FROM COMMAND LINE
beta = int(sys.argv[1])
tfinal = float(sys.argv[2])
n = int(sys.argv[3])
chi = float(sys.argv[4])
gamma = float(sys.argv[5])
logstep = float(sys.argv[6])

color = ['b','r', 'g', 'k']    
cortype = ['cor2p', 'cor2r','cor2rp', 'cor2pr']
l = 0
for cor in cortype:
    if cor == 'cor2p':
        cor_name = 'S_{22}'
        
    if cor == 'cor2r':
        cor_name = 'S_{11}'
        
    if cor == 'cor2pr':
        cor_name = 'S_{21}'
        
    if cor == 'cor2rp':
        cor_name = 'S_{12}'
    l = l + 1
    basename = 'FPU_%s_n_%d_beta_%d_chi_%f_gamma_%f' %(cor,n,beta,chi,gamma) # that is the same name that we gave to the file 
    for f in glob.glob('Tot_*%s*.dat' %(basename)):                          # produced by FPUT_short_range.py
        data = np.loadtxt(f)
        #Time scale
        timer = spaziatura_log(tfinal,logstep)
        t = len(timer)
        
        
        # Plot of the correlation functions
        color_counter = -1
        for k in np.arange(0,t, step = int(t/4)):
        
            color_counter = color_counter + 1
            xfine = np.linspace(-0.5,0.5,n)
            plt.plot(xfine,data[k,:],alpha = 0.8, label = r't = %0.f' %timer[k], color = color[color_counter])
                    
        plt.xlabel(r'$j/N$')
        plt.ylabel(r'$%s$' %cor_name)
        plt.title(r'$%s$, $\chi=%0.2f$, $\gamma =%0.3f$' %(cor_name,chi,gamma))
        plt.legend(loc = 1)
            
        plt.savefig('Simple_cor_%s_chi_%0.3f_gamma_%0.3f.png' %(ex_number,cor,chi,gamma))
        plt.close()
        
        # Plot of the scaling of the pick at j =0, if it is not present comment
        
        maximum_peak = np.amax(np.abs(data[:, -100 + int(n*0.5): 100 + int(n*0.5)]), axis=1)

        fit_param,error = opt.curve_fit(func,timer[10:] , maximum_peak[10:])
        plt.loglog(timer[1:],maximum_peak[1:], 'k.', label = 'Central Peak', markersize=30)
        plt.loglog(timer[1:], func(timer[1:], fit_param[0],fit_param[1]), 'r-',label = r'Interpolation: $%.3f t^{%.3f}$' %(fit_param[1],fit_param[0]),  linewidth = 4)
        plt.xlabel(r'$t$')
        plt.ylabel('|Central Peak $%s$|' %cor_name)
        plt.title(r'Scaling central peak $%s$, $\chi=%0.2f$, $\gamma =%0.3f$' %(cor_name,chi,gamma))
        plt.legend(loc = 1)
        plt.xlim(100,1000)
        plt.savefig('Interpolation_central_peak_%s_chi_%0.3f_gamma_%0.3f_logplot.png' %(ex_number,cor,chi,gamma))
        plt.close()
        
        
        
        # Plot of the extreme peak assuming there are no other peak in the correlation, except the central one
        round_value = 300
                
        decay = 2/3 # guess of the decay 1/3 for weak nonlinearity, 2/3 for strong non linearity
        
        xfine = np.linspace(-timer[-1]**(decay),timer[-1]**(decay),num = round_value)
        color_counter = -1
        for k in np.arange(0,t,step = int(t/4)):
            color_counter = color_counter + 1 

            max_pos = np.argmax(data[:,int(n/2) + 100:], axis=1) + int(n/2) + 100

                    
            effective_position_minus = int(max_pos[k] - 0.5*(round_value))
            effective_position_plus = int(max_pos[k] + 0.5*(round_value))
                    
            plt.plot(xfine/(timer[k]**(decay)),(timer[k]**decay)*data[k,effective_position_minus:effective_position_plus],linewidth = 3, alpha = 0.8, label = r't = %0.f' %timer[k], color = color[color_counter])
                    
            plt.xlabel(r'$(j-\xi_0t)/(Nt^\frac{2}{3})$')
            plt.ylabel(r'$t^\frac{2}{3}%s$' %cor_name)
            plt.title(r'Scaling extreme peak $%s$, $\chi=%0.2f$, $\gamma =%0.3f$' %(cor_name,chi,gamma))
            plt.legend(loc = 1)
            plt.xlim(-1,1)
            plt.savefig('Scaling_extreme_peak_%s_chi_%0.3f_gamma_%0.3f.png' %(ex_number,cor,chi,gamma))
            plt.close()
