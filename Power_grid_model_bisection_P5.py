# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 09:10:33 2022

@author: Paul

Script for calculating critical boundaries of power grid model using 
bisection method 
"""

#import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.io import savemat
from scipy.io import loadmat


def BusModel_4D_PG_IVP(t,x,Q_shift,r):

    """
    3-bus model for a power grid
    """
    
  
    # Variables
    delta_m = x[0]
    w_m = x[1]
    delta = x[2]
    V = x[3]
    Q_1 = x[4]
    
    # Generator parameters
    M = 0.3
    d_m = 0.05
    P_m = 1
    
    # Load parameters
    K_pw = 0.4
    K_pv = 0.3
    K_qw = -0.03
    K_qv = -2.8
    K_qv2 = 2.1
    T = 8.5
    P_0 = 0.6
    Q_0 = 1.3
    
    # Network parameters
    V_m = 1
    Y_m = 5
    theta_m = -5*np.pi/180
    V_0prime = 2.5
    Y_0prime = 8
    theta_0prime = -12*np.pi/180
    
    P_l = -V_0prime*V*Y_0prime*np.sin(delta + theta_0prime) - V_m*V*Y_m*np.sin(delta - delta_m + theta_m) + (Y_0prime*np.sin(theta_0prime) + Y_m*np.sin(theta_m))*(V**2)
    Q_l = -(-V_0prime*V*Y_0prime*np.cos(delta + theta_0prime) - V_m*V*Y_m*np.cos(delta - delta_m + theta_m) + (Y_0prime*np.cos(theta_0prime) + Y_m*np.cos(theta_m))*(V**2))
    P_e = -V_m*V*Y_m*np.sin(delta - delta_m - theta_m) - (V_m**2)*Y_m*np.sin(theta_m)
    
    z = np.zeros(5)

    z[0] = w_m
    z[1] = (-d_m*w_m + P_m - P_e)/M
    z[2] = (-K_qv*V - K_qv2*(V**2) + Q_l - Q_0 - Q_1)/K_qw
    z[3] = (K_pw*K_qv2*(V**2) + (K_pw*K_qv - K_qw*K_pv)*V + K_pw*(Q_0 + Q_1 - Q_l) - K_qw*(P_0 + P_1 - P_l))/(T*K_qw*K_pv)
    z[4] = -(Q_shift*r*np.tanh(r*t)/np.cosh(r*t))#*(t<0)+0*(t>=0)#0
    
    return(z)
    
def event(t,x,Q_shift,r):
    return x[3]
    
event.terminal = True

## Initial conditions

delta_m = 0.1119
w_m = 0
delta = -0.097#-0.085
V = 1.5149
Q_1 = 5
P_1 = 5

# Q_shifts = np.linspace(5.9,6,11)
rs = np.geomspace(1,1E3,101)#(8E0,1E2,41)#(1E0,1E1,33)#(1E0,1E3,101)#,101)

## Set up for initial value problem solver

tspan = [-5,5]#[-50,2000]
h = 0.001
t = np.arange(tspan[0],tspan[1],h)

Nmax = 12#20
Q_shifts = np.zeros(len(rs))

for i in range(len(rs)):
    
    Q_shifta = 2
    Q_shiftb = 13
    
    N = 0
    
    while N<Nmax:
        
        Q_shiftc = (Q_shiftb + Q_shifta)/2
        
        sol = solve_ivp(BusModel_4D_PG_IVP,tspan,[delta_m, w_m, delta, V, Q_1],args=(Q_shiftc,rs[i]),events=event,max_step=h)
        
        # if sol.y[3,-1]>1E-6 and sol.y[2,-1]<5:
        if sol.y[3,-1]>1E-6:
            Q_shifta = Q_shiftc
        else:
            Q_shiftb = Q_shiftc
        
        N = N+1

    Q_shifts[i] = Q_shifta
    print (i)

plt.figure()
plt.plot(Q_shifts,rs)

# Create dictionary of data to save                
mdic = {'Q_1': Q_1, 'Q_shifts': Q_shifts, 'rs': rs, 'Nmax': Nmax, 'Q_shifta': Q_shifta, 'Q_shiftb': Q_shiftb, 'tspan': tspan, 'delta_m': delta_m, 'w_m': w_m, 'delta': delta, 'V': V}

############# Save data as a matifle - Specify filename #############
savemat('Power_grid_model_monotone_forcing_tipping_fast_rates_data.mat', mdic)
