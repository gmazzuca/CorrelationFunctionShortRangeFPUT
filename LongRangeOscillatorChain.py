''' Thi script compute the correlation functions for several example of the short range oscillator chain, i.e. for a family
of Hamiltonian systems of the form H = \sum_{j=0}^{N-1} (p_j^2/2 + \sum_{s=1}^m k_s(q_{j+s} - q_j)^2/2) = p^Tp/2 + q^TAq/2, where m<<N, we call k_1,..,k_m
the intensity of the springs and A is a circulant matrix.
'''

#################### Imports ##############
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.pylab as pylab

################## parameters for plots ################
params = {'legend.fontsize': 30,
          'figure.figsize': (17, 13),
         'axes.labelsize': 45,
         'axes.titlesize':50,
         'xtick.labelsize':25,
          'ytick.labelsize':25}
pylab.rcParams.update(params)



###################### external functions ####################
def pearcey(s):
    '''This function returns Pearcey integral absolute value squared, this functions is taken from [1] '''
    rate = 0.01
    x = np.linspace(-4.0, 4.0, 80)
    
    y = rate*(4*x**3+s)
    z = x+1j*y

    f = z**4+s*z
    g = np.exp(1j*f)

    # ∫f(z)dz = ∫f(z)(dz/dx)dx
    dz = 1.0+1j*rate*(12*x**2)
    I = integrate.simps(g*dz, x)

    return np.abs(I)**2


###################################################################
###################################################################
######################### MAIN ####################################
###################################################################
###################################################################

eps = 0.000001 # to avoid round off error
particles = 10001 # size of the chain
beta = 1 # inverse of temperature
timesample = 3 # number of snaps
color = ['b','r', 'g'] # colors for the plot
betainv = 1/beta # temperature
timer = np.linspace(500,1500,timesample) # time range
jfine = np.arange(particles)- int(particles*0.5) # particles range
xfine = np.linspace(0+eps,1-eps,particles) # spectrum



# Several examples of Short Range Harmonic Chain
# k1,...,km are the intensity of the springs
# xi is the ratio j/t where j is the position of the particles and t is the time
# v0,v1 are the velocities of peaks, v0 is the fastest one
# lambda0 is the Airy parameter
# lambda1 is the Pearcey parameter
# h is the dispersive relation
# t1,..,tm are the entries of T, the local square root of A
# theta is a function that describe the arguments of the eigenvalues of the matrix T

for l in range(3):
            
            
    if l == 0:
        # classical harmonic oscillator
        k1 = 1
        xi = 1;
        decay = 1/4
        v0= np.sqrt(k1)
        lambda0 = 0.5
        h = lambda x : 2*k1*(1-np.cos(2*np.pi*x) )
        h_prime = lambda x : 2*2*np.pi*k1*(np.sin(2*np.pi*x))
        theta = lambda x : np.arctan(1/np.tan(np.pi*x))
        lamda1 = 0
        xlimit = 2100
            
    if l == 1:
        ## example with a traveling peak with decay t^{-1/4} and a traveling peak with decay t^{-1/3}
        k1 = 1
        k2 = 1/8
        k3 = 7/72
        xi = 0.5/np.sqrt(2)
        decay = 1/4
        v0 = np.sqrt(k1 + 4*k2+ 9*k3)
        lambda0 = 0.959038018341891
        h = lambda x: 2*k1*(1-np.cos(2*np.pi*x)) + 2*k2*(1-np.cos(4*np.pi*x)) + 2*k3*(1-np.cos(6*np.pi*x))
        h_prime = lambda x : 2*2*np.pi*k1*(np.sin(2*np.pi*x)) + 2*4*np.pi*k2*(np.sin(4*np.pi*x)) + 2*6*np.pi*k3*(np.sin(6*np.pi*x))
        lambda1 =   0.781648328813151
        xlimit = 2800
        t0 = 1.21422
        t1 = -0.967414;
        t2 = -0.16674;
        t3 = -0.0800694;
        theta = lambda x : np.arctan(-(t1*np.sin(2*np.pi*x) + t2*np.sin(4*np.pi*x) + t3*np.sin(6*np.pi*x))/(t0 +t1*(np.cos(2*np.pi*x)) + t2*(np.cos(4*np.pi*x)) + t3*(np.cos(6*np.pi*x))))
        
    elif l == 2:
        ## example with a fixed peak with decay t^{-1/4} and a traveling peak with decay t^{-1/3}
        k1 = 1
        k2 = 0.25
        xi = 0
        decay = 1/4
        v0 = np.sqrt(k1 + 4*k2)
        lambda0 = 0.761708
        h = lambda x : 2*k1*(1-np.cos(2*np.pi*x)) + 2*k2*(1-np.cos(4*np.pi*x))
        h_prime = lambda x : 2*2*np.pi*k1*(np.sin(2*np.pi*x)) + 2*4*np.pi*k2*(np.sin(4*np.pi*x)) 
        lambda1 = 0.5
        xlimit = 2800
        t0 = 0.5*( + np.sqrt(k1) +np.sqrt(k1+4*k2))
        t1 = -np.sqrt(k1)
        t2 = -0.5*( - np.sqrt(k1) + np.sqrt(k1+4*k2))
        theta = lambda x : np.arctan(-(t1*np.sin(2*np.pi*x) + t2*np.sin(4*np.pi*x))/(t0 + t1*np.cos(2*np.pi*x) + t2*np.cos(4*np.pi*x)))
        
        
    
       
       
    # Plot of the dispersive relations
    plt.plot(xfine, np.sqrt(h(xfine)), 'b', linewidth = 3)
    plt.title('Dispersive relation')
    plt.xlabel(r'$k$')
    plt.ylabel(r'$f$', rotation = 0)
    plt.savefig('dispersive_ex_%d.png' %l)
    plt.close()
    
    # Plot of the derivative of the dispersive relations
    
    plt.plot(xfine, 0.5*h_prime(xfine)/np.sqrt(h(xfine)), 'b', linewidth = 3)
    plt.title('Derivative of the dispersive relation')
    plt.xlabel(r'$k$')
    plt.ylabel(r'$f^{\prime}$', rotation = 0)
    plt.savefig('Derivative_dispersive_ex_%d.png' %l)
    plt.close()

    # Compute the correlations 

    S22 = lambda j,t :  np.trapz(  np.cos(np.sqrt(h(xfine))*t)*np.cos(2*np.pi*xfine*j), x = xfine)
    S21 = lambda j,t :  np.trapz(  -np.sin(np.sqrt(h(xfine))*t)*np.cos(2*np.pi*xfine*j - theta(xfine)), x = xfine)
    S12 = lambda j,t :  np.trapz(  np.sin(np.sqrt(h(xfine))*t)*np.cos(2*np.pi*xfine*j + theta(xfine)), x = xfine)


    # Correlation finals 
    S22plot = np.zeros((timesample,particles))
    S21plot = np.zeros((timesample,particles))
    S33plot = np.zeros((timesample,particles))
    S12plot = np.zeros((timesample,particles))

    for i in np.arange(particles):
        for t in np.arange(timesample):
            S22plot[t,i] = S22(jfine[i],timer[t])
            S21plot[t,i] = S21(jfine[i],timer[t])
            S12plot[t,i] = S12(jfine[i],timer[t])
            S33plot[t,i] = (S21plot[t,i]*S21plot[t,i] + S12plot[t,i]*S12plot[t,i])*0.5 + S22plot[t,i]*S22plot[t,i]
            
    
    fname = 'example_%d.dat' %l
    # save the data for the correlations
    np.savetxt('S22_%s' %fname, S22plot)
    np.savetxt('S21_%s' %fname, S21plot)
    np.savetxt('S33_%s' %fname, S33plot)
    
    
    # Plot just the correlation functions
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()
    
    
    for k in np.arange(timesample):

        ax1.plot(jfine, betainv*S22plot[k,:], color[k], alpha = 0.8, label ='t = %.0f'%timer[k],  linewidth = 2 )
        ax1.set_title(r'Correlation $S_{11},S_{22}$')
        ax1.set_ylabel(r'$S_{11},S_{22}$')
        ax1.set_xlabel(r'$j$')
        ax1.set_xlim(-xlimit,xlimit)
        
        
        ax2.plot(jfine, betainv*S21plot[k,:], color[k], alpha = 0.8,label ='t = %.0f'%timer[k],  linewidth = 2)
        ax2.set_title(r'Correlation $S_{21}$')
        ax2.set_ylabel(r'$S_{21}$')
        ax2.set_xlabel(r'$j$')
        ax2.set_xlim(-xlimit,xlimit)
        
        
        ax3.plot(jfine, betainv*S33plot[k,:], color[k], alpha = 0.8,label ='t = %.0f'%timer[k],  linewidth = 2)
        ax3.set_title(r'Correlation $S_{33}$')
        ax3.set_ylabel(r'$S_{33}$')
        ax3.set_xlabel(r'$j$')
        ax3.set_xlim(-xlimit,xlimit)

    ax1.legend(loc =1 )
    ax2.legend(loc =1 )
    ax3.legend(loc =1 )
    fig1.savefig('S11_Ex_%d_simp_cor.png' %l)
    fig2.savefig('S21_Ex_%d_simp_cor.png' %l)
    fig3.savefig('S33_Ex_%d_simp_cor.png' %l)
    plt.close('all')

    if l > 0:
        # Scaling Pearcey Peak
        lower = 400
        xlaterale = np.linspace(-lower,lower,2*lower)
        fig1, ax1 = plt.subplots()
        fig2, ax2 = plt.subplots()
        fig3, ax3 = plt.subplots()
        
        
        
        for k in np.arange(timesample):
            xmin = int(particles*0.5) + int(xi*timer[k]) - lower
            xmax = int(particles*0.5) + int(xi*timer[k]) + lower
            xlimit = 30
            ax1.plot(xlaterale/(lambda1*np.power(timer[k],1/4)), np.power(timer[k],decay)*betainv*S22plot[k,xmin:xmax], color = color[k], linewidth = 3, alpha = 0.8,label ='t = %.0f'%timer[k])
            ax1.set_title(r'Pearcey scaling $S_{11}$')
            ax1.set_ylabel(r'$t^{\frac{1}{4}}S_{11}$')
            ax1.set_xlabel(r'$(j-v^*t)/(\lambda^*t^\frac{1}{4})$')
            ax1.set_xlim(-xlimit,xlimit)
            
            ax2.plot(xlaterale/(lambda1*np.power(timer[k],1/4)), np.power(timer[k],decay)*betainv*S21plot[k,xmin:xmax],  linewidth = 3, color = color[k], alpha = 0.8,label ='t = %.0f'%timer[k])
            ax2.set_title(r'Pearcey scaling  $S_{21}$')
            ax2.set_ylabel(r'$t^{\frac{1}{4}}S_{12}$')
            ax2.set_xlabel(r'$(j-v^*t)/(\lambda^*t^\frac{1}{4})$')
            ax2.set_xlim(-xlimit,xlimit)
            
            ax3.plot(xlaterale/(lambda1*np.power(timer[k],1/4)), np.power(timer[k],2*decay)*betainv*S33plot[k,xmin:xmax], linewidth = 3,  color =color[k], alpha = 0.8,label ='t = %.0f'%timer[k])
            ax3.set_title(r'Pearcey scaling  $S_{33}$')
            ax3.set_ylabel(r'$t^{\frac{1}{2}}S_{33}$')
            ax3.set_xlabel(r'$(j-v^*t)/(\lambda^*t^\frac{1}{4})$')
            ax3.set_xlim(-xlimit,xlimit)
            
        if (l == 2):
            # Plot Pearcy integral
            p = np.zeros(len(xlaterale))
            ifine = np.arange(-lower, lower)
            for ll,s in enumerate(ifine):
                p[ll] = pearcey(s/(lambda1*np.power(timer[k],1/4)))
            ax3.plot(xlaterale/(lambda1*np.power(timer[k],1/4)), p*(betainv**2)/(4*np.pi*np.pi*lambda1*lambda1), linewidth = 3, color ='k',label ='Pearcey')
            ax3.set_xlim(-xlimit,xlimit)


        ax1.legend(loc =1 )
        ax2.legend(loc =1 )
        ax3.legend(loc =1 )
        fig1.savefig('S11_Ex_%d_scaling_center_pick.png' %l)
        fig2.savefig('S21_Ex_%d_scaling_center_pick.png' %l)
        fig3.savefig('S33_Ex_%d_scaling_center_pick.png' %l)
        plt.close('all')


    # Scaling Airy Peak
    
    
    lower = 500
    xlaterale = np.linspace(-timer[2]**(1/3)*lambda0,timer[2]**(1/3)*lambda0,2*lower)
    special_airy = lambda x : betainv*sc.special.airy(x)[0]/(2*lambda0)
    fig1, ax1 = plt.subplots()
    fig2, ax2 = plt.subplots()
    fig3, ax3 = plt.subplots()

    for k in np.arange(timesample):
        xmin = int(particles*0.5) + int(v0*timer[k]) - lower
        xmax = int(particles*0.5) + int(v0*timer[k]) + lower
        
        
        ax1.plot(xlaterale/(lambda0*np.power(timer[k],1/3)), np.power(timer[k],1/3)*betainv*S22plot[k,xmin:xmax],linewidth = 3, color = color[k], alpha = 0.8,label ='t = %.0f'%timer[k])
        ax1.set_title(r'Airy scaling  $S_{11}$')
        ax1.set_ylabel(r'$t^{\frac{1}{3}}S_{11}$')
        ax1.set_xlabel(r'$(j-v_0t)/(\lambda_0t^{\frac{1}{3}})$')
        ax1.set_xlim(-0.25,0.25)
        
        
        ax2.plot(xlaterale/(lambda0*np.power(timer[k],1/3)), np.power(timer[k],1/3)*betainv*S21plot[k,xmin:xmax],linewidth = 3,color = color[k], alpha = 0.8,label ='t = %.0f'%timer[k])
        ax2.set_title(r'Airy scaling  $S_{21}$')
        ax2.set_ylabel(r'$t^{\frac{1}{3}}S_{21}$')
        ax2.set_xlabel(r'$(j-v_0t)/(\lambda_0t^{\frac{1}{3}})$')
        ax2.set_xlim(-0.25,0.25)
        
        ax3.plot(xlaterale/(lambda0*np.power(timer[k],1/3)), np.power(timer[k],2/3)*betainv*S33plot[k,xmin:xmax],linewidth = 3, color =color[k], alpha = 0.8,label ='t = %.0f'%timer[k])
        ax3.set_title(r'Airy scaling  $S_{33}$')
        ax3.set_ylabel(r'$t^{\frac{2}{3}}S_{33}$')
        ax3.set_xlabel(r'$(j-v_0t)/(\lambda_0t^{\frac{1}{3}})$')
        ax3.set_xlim(-0.25,0.25)


    
    ifine = np.arange(-lower, lower)
    ax1.plot(xlaterale/(lambda0*np.power(timer[2],1/3)), special_airy(ifine/(lambda0*timer[2]**(1/3))),linewidth = 3, color ='k', label = 'Airy') 
    ax2.plot(xlaterale/(lambda0*np.power(timer[2],1/3)), -special_airy(ifine/(lambda0*timer[2]**(1/3))),linewidth = 3,  color ='k', label = 'Airy') 
    ax3.plot(xlaterale/(lambda0*np.power(timer[2],1/3)), 2*special_airy(ifine/(lambda0*timer[2]**(1/3)))*special_airy(ifine/(lambda0*timer[2]**(1/3))),linewidth = 3,  color ='k', label = 'Airy')
    ax1.legend(loc =1 )
    ax2.legend(loc =1 )
    ax3.legend(loc =1 )
    fig1.savefig('S11_Ex_%d_scaling_extreme_pick.png' %l)
    fig2.savefig('S21_Ex_%d_scaling_extreme_pick.png' %l)
    fig3.savefig('S33_Ex_%d_scaling_extreme_pick.png' %l)
    plt.close('all')


'''
[1] Dan Piponi, "Pearcey Integral", GitHub Repository  URL: https://gist.github.com/dpiponi/9176c7f6bf32803e9b2bf6e8c0b93ab5
'''