# -*- coding: utf-8 -*-
"""
Created on Thu May  4 10:43:00 2023

@author: Paul

Figure 3: Tipping for a nonlinear shift ramp and return profiles in the tilted fold example
"""

#import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy import interpolate
from scipy.integrate import solve_ivp
from scipy.io import savemat
from scipy.io import loadmat
from matplotlib import rc
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
import seaborn as sns

plt.rcParams["font.family"] = "serif"
fontsize = 12
labelsize=10



def f(x, p, A, s):
    """
    Drift
    """  
    f = -x*( ( x - A - s*p )**2  + p)   
    return (f)

def Forcing(t, p0, Deltap, r):
    """
    External forcing ramp
    """      
    p = (p0 + Deltap/np.cosh(r*t))*(t<=0) + (p0+Deltap)*(t>0) 
    return (p)

def Forcing2(t, p0, Deltap, r):
    """
    External forcing return
    """      
    p = (p0 + Deltap/np.cosh(r*t))
    return (p)

## Time parameters
tstart = -500                              # Start time
tend = 10                              # End time
dt = 0.001                                # Spacing between time intervals
n = int((tend - tstart)/dt)             # Number of time intervals
t = np.linspace(tstart, tend, n+1)      # Time values

## System parameters
s = 4                               # Tilt parameter
A = 3.2                             # Separation parameter


## Forcing parameters
# Forcing rate parameters
r1      = 0.5
r2      = 0.9
r3      = 1.7
# Max shift distance
Deltap = 0.48
# Start value of forcing
p0 = -0.5

p1 = Forcing(t, p0, Deltap, r1)
p2 = Forcing(t, p0, Deltap, r2)
p3 = Forcing(t, p0, Deltap, r3)
p4 = Forcing2(t, p0, Deltap, r1)
p5 = Forcing2(t, p0, Deltap, r2)
p6 = Forcing2(t, p0, Deltap, r3)


## Initialise variable 
x = np.zeros(n+1)
x2 = np.zeros(n+1)
x3 = np.zeros(n+1)
x4 = np.zeros(n+1)
x5 = np.zeros(n+1)
x6 = np.zeros(n+1)
x[0] = A + s*p0 + np.sqrt(-p0)
x2[0] = A + s*p0 + np.sqrt(-p0)
x3[0] = A + s*p0 + np.sqrt(-p0)
x4[0] = A + s*p0 + np.sqrt(-p0)
x5[0] = A + s*p0 + np.sqrt(-p0)
x6[0] = A + s*p0 + np.sqrt(-p0)


for i in range(n):
    ## Solve ODE with Euler method
    x[i+1] = x[i] + dt*f(x[i], p1[i], A, s)
    x2[i+1] = x2[i] + dt*f(x2[i], p2[i], A, s)
    x3[i+1] = x3[i] + dt*f(x3[i], p3[i], A, s)
    x4[i+1] = x4[i] + dt*f(x4[i], p4[i], A, s)
    x5[i+1] = x5[i] + dt*f(x5[i], p5[i], A, s)
    x6[i+1] = x6[i] + dt*f(x6[i], p6[i], A, s)


## Initialise figure
fig=plt.figure(figsize=(11,4.5))
gs=GridSpec(3,2,hspace=0,height_ratios=[1,0.05,1])#,0.07,1.8])

ax1=fig.add_subplot(gs[0,0]) 
ax2=fig.add_subplot(gs[2,0],sharex=ax1)
ax3=fig.add_subplot(gs[:,1])

## Plotting
ax1.plot(t[t<0],p1[t<0],c=[25/255,101/255,176/255],ls='--')
ax1.plot(t[t<0],p2[t<0],c=[238/255,128/255,38/255],ls='--')
ax1.plot(t[t<0],p3[t<0],c=[174/255,118/255,163/255],ls='--')
ax1.plot(t[t>=0],p1[t>=0],c='k',ls='--')
ax1.plot(t,p4,c=[25/255,101/255,176/255])
ax1.plot(t,p5,c='tab:orange')
ax1.plot(t,p6,c='tab:purple')
ax1.set_ylabel('External forcing',fontsize=fontsize,labelpad=10)
ax1.set_xticklabels([])
ax1.set_yticklabels([])

ax2.plot(t,x,c=[25/255,101/255,176/255],ls='--')
ax2.plot(t,x2,c=[238/255,128/255,38/255],ls='--')
ax2.plot(t,x3,c=[174/255,118/255,163/255],ls='--')
ax2.plot(t,x4,c=[25/255,101/255,176/255])
ax2.plot(t,x5,c=[238/255,128/255,38/255])
ax2.plot(t,x6,c=[174/255,118/255,163/255])
ax2.set_xlabel('Time',fontsize=fontsize,labelpad=10)
ax2.set_ylabel('System state',fontsize=fontsize,labelpad=10)
ax2.set_xticklabels([])
ax2.set_yticklabels([])
ax2.set_xlim(-10,10)





## Import data for tipping diagram
mat1 = loadmat('Tilted_SN_critical_boundary_mon_forcing_data_p05_slow_rates.mat')
mat2 = loadmat('Tilted_SN_critical_boundary_mon_forcing_data_p05_fast_rates.mat')
mat3 = loadmat('Tilted_SN_critical_boundary_nonmon_forcing_data_p05.mat')

Deltap_shifts = np.concatenate((mat1['Deltaps'][0,:],mat2['Deltaps'][0,1:]))
rs = np.concatenate((mat1['rs'][0,:],mat2['rs'][0,1:]))
Deltap_returns = mat3['Deltaps'][0,:]
rs_returns = mat3['rs'][0,:]

f = interpolate.interp1d(rs, Deltap_shifts)
f2 = interpolate.interp1d(rs_returns, Deltap_returns)

x = np.geomspace(0.01,1000,10001)
x2 = np.geomspace(0.01,10,10001)
x3 = np.geomspace(10,1000,10001)
y1 = f(x)
y2 = f2(x2)
y3 = f(x2)
y4 = f(x3)


## Plotting critical boundaries

ax3.plot(y1,x,c='k',linestyle='dashed')
ax3.plot(y2,x2,c='k')
ax3.plot([0.5,0.5],[0.01,1000],'k',lw=0.5)
BI = 0.29091797
ax3.plot([BI,BI],[0.01,1000],'k',lw=0.5)

ax3.fill_betweenx(x2,y2,y3,color=[78/255,178/255,101/255],edgecolor='none',alpha=0.4)
ax3.fill_betweenx(x3,y4,1,color=[78/255,178/255,101/255],edgecolor='none',alpha=0.4)
ax3.fill_betweenx(x2,y2,1,color=[220/255,5/255,12/255],edgecolor='none',alpha=0.4)
ax3.set_xlim(0.25,0.6)
ax3.set_ylim(0.1,50)
ax3.set_yscale('log')
ax3.set_xlabel('Peak change in external forcing',fontsize=fontsize,labelpad=10)
ax3.set_ylabel('Rate parameter of external forcing',fontsize=fontsize,labelpad=10)
ax3.set_xticklabels([])
ax3.set_yticklabels([])


sns.despine()
ax3.spines["top"].set_visible(True)
ax3.spines["right"].set_visible(True)

ax3.plot(Deltap,r1,Marker='o',ms=10,c=[25/255,101/255,176/255])
ax3.plot(Deltap,r2,Marker='o',ms=10,c=[238/255,128/255,38/255])
ax3.plot(Deltap,r3,Marker='o',ms=10,c=[174/255,118/255,163/255])
ax3.text(0.492, 30, "Fold", size=labelsize, rotation=90, ha="center", va="center")
ax3.text(0.285, 1, "Basin instability boundary", size=labelsize, rotation=90, ha="center", va="center")

ax1.text(0.03, 0.86, 'a)', transform=ax1.transAxes,size=12,fontsize=fontsize, weight='bold')
ax2.text(0.03, 0.86, 'b)', transform=ax2.transAxes,size=12,fontsize=fontsize, weight='bold')
ax3.text(0.03, 0.93, 'c)', transform=ax3.transAxes,size=12,fontsize=fontsize, weight='bold')

fig.tight_layout()