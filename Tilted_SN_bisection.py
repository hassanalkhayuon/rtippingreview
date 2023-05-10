# -*- coding: utf-8 -*-
"""
Created on Fri Feb  4 15:20:13 2022

@author: pdlr201

Bisection method to calculate critical boundaries for tipping diagram Fig. 3
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import seaborn as sns
from matplotlib.gridspec import GridSpec
import matplotlib.patches as patches

def f(x, p, A, s):
    """
    Drift
    """  
    f = -x*( ( x - A - s*p )**2  + p)   
    return (f)

def Forcing(t, p0, Deltap, r):
    """
    External forcing
    """      
    p = (p0+Deltap/np.cosh(r*t))#*(t<0) + (p0+Deltap)*(t>=0)
    return (p)

# System parameters
s = 4
A = 3                             # Inverse timescale parameter
kappa = 3                               # Curvature parameter


p0 = -0.4                                 # Start value of forcing

rs = np.geomspace(1E-1,1E3,101)#0#np.geomspace(1E-4,1E-0,101)
err = 0.001#0.0001
Deltaps = np.zeros(len(rs))


for j in range(len(rs)):
    
    Deltapa = 0.1#0.02#0
    Deltapb = 0.6#0.06#1.6
    
    Deltapc = (Deltapb + Deltapa)/2
    
    eps = np.abs(Deltapb - Deltapa)
    
    N = 0

    tspan = [-8/rs[j],np.maximum(40/rs[j],8)]
    dt = 0.0001#0.000005/epsilonrs[j]
    n = int((tspan[1]-tspan[0])/dt)
    t = np.linspace(tspan[0], tspan[1], n)
    
    while eps>err:
        # Initialise variable 
        x = np.zeros(n+1)
        
        x[0] = A + s*p0 + np.sqrt(-p0)

        for i in range(n):
            # Solve ODE with Euler method
            x[i+1] = x[i] + dt*f(x[i], Forcing(t[i], p0, Deltapc, rs[j]), A, s)
            
        if x[-1] < 0.5:
            Deltapb = Deltapc
        else:
            Deltapa = Deltapc
            
        Deltapc = (Deltapb + Deltapa)/2
        
        eps = np.abs(Deltapb - Deltapa)
        N = N+1
    
    Deltaps[j] = Deltapc
    print (j, N)


plt.figure()
plt.plot(Deltaps,rs)

