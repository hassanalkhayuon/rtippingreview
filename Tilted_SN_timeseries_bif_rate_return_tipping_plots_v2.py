# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 16:52:49 2023

@author: pdlr201

ESD Fig 2:  Illustration of bifurcation-induced, rate-induced, and return tipping
"""

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import seaborn as sns
from matplotlib.gridspec import GridSpec
import matplotlib.patches as patches

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
    External forcing
    """    
    # p = p0 + Deltap*(np.tanh(r*t) + 1)/2   
    p = (p0 + Deltap/np.cosh(r*t))*(t<=0) + (p0+Deltap)*(t>0) 
    return (p)

def Forcing2(t, p0, Deltap, r, r2):
    """
    External forcing
    """    
    # p = p0 + Deltap*(np.tanh(r*t) + 1)/2   
    p = (p0 + Deltap/np.cosh(r*t))*(t<=0) + (p0 + Deltap/np.cosh(r2*t))*(t>0) 
    return (p)

## Time parameters
tstart = -500                              # Start time
tend = 100                              # End time
dt = 0.01                                # Spacing between time intervals
n = int((tend - tstart)/dt)             # Number of time intervals
t = np.linspace(tstart, tend, n+1)      # Time values

## System parameters
# Tilt parameters
s1 = 4
s2 = -4
#Separation parameters
A1 = 3.2   
A2 = 1.5                           


## Forcing parameters
# Forcing rate parameters
r1      = 0.05
r2      = 1
r3      = 3
# Max shift distance
Deltap1 = 0.505
Deltap2 = 0.4
# Start value of forcing
p0 = -0.5                                 


p1 = Forcing(t, p0, Deltap1, r1)
p2 = Forcing(t, p0, Deltap2, r2)
p3 = Forcing(t, p0, Deltap2, r3)
p4 = Forcing(t, p0+Deltap2, -Deltap2, r2)
p5 = Forcing(t, p0+Deltap2, -Deltap2, r3)


## Equilibria
pp = np.linspace(-0.9, 0.2, 601)
xeqplus1 = A1 + s1*pp + np.sqrt(-pp)
xeqminus1 = A1 + s1*pp - np.sqrt(-pp)
xeqplus2 = A2 + s2*pp + np.sqrt(-pp)
xeqminus2 = A2 + s2*pp - np.sqrt(-pp)

idx1 = np.nanargmin(np.abs(xeqminus1-(A1 + s1*p0 + np.sqrt(-p0))))
idx2 = np.nanargmin(np.abs(xeqminus2-(A2 + s2*(p0+Deltap2) + np.sqrt(-(p0+Deltap2)))))
BI1 = pp[idx1]
BI2 = pp[idx2]

## Initialise variable 
x = np.zeros(n+1)
x2 = np.zeros(n+1)
x3 = np.zeros(n+1)
x4 = np.zeros(n+1)
x5 = np.zeros(n+1)
x[0] = A1 + s1*p0 + np.sqrt(-p0)
x2[0] = A1 + s1*p0 + np.sqrt(-p0)
x3[0] = A1 + s1*p0 + np.sqrt(-p0)
x4[0] = A2 + s2*(p0+Deltap2) + np.sqrt(-(p0+Deltap2))
x5[0] = A2 + s2*(p0+Deltap2) + np.sqrt(-(p0+Deltap2))


## Initialise figure
fig=plt.figure(figsize=(13.5,5.25))

gs=GridSpec(2,3,height_ratios=[1,2])

ax1=fig.add_subplot(gs[0,0]) 
ax2=fig.add_subplot(gs[0,1],sharey=ax1)
ax3=fig.add_subplot(gs[0,2],sharex=ax2,sharey=ax1)
ax4=fig.add_subplot(gs[1,0])
ax5=fig.add_subplot(gs[1,1],sharex=ax4,sharey=ax4)
ax6=fig.add_subplot(gs[1,2],sharex=ax4)

ax1.set_xlim(-150, 50)
ax1.set_ylim(-0.55, 0.1)
ax2.set_xlim(-6, 2)
ax1.set_xlabel('Time', fontsize = fontsize)
ax1.set_ylabel('External forcing', fontsize = fontsize)
ax2.set_xlabel('Time', fontsize = fontsize)
ax3.set_xlabel('Time', fontsize = fontsize)
ax1.tick_params(axis='both', which='major', labelsize=10)
ax2.tick_params(axis='both', which='major', labelsize=10)
ax3.tick_params(axis='both', which='major', labelsize=10)
ax1.set_xticks([])
ax1.set_yticks([])
ax2.set_xticks([])
ax2.set_yticks([])
ax3.set_yticks([])

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['top'].set_visible(False)

## Plotting
ax1.plot([-150, 50], [0, 0], 'k--', lw=1)
ax2.plot([-150, 50], [0, 0], 'k--', lw=1)
ax3.plot([-150, 50], [0, 0], 'k--', lw=1)

ax1.text(0.05, 0.88, 'Fold', transform=ax1.transAxes,size=12)
ax2.text(0.05, 0.88, 'Fold', transform=ax2.transAxes,size=12)
ax3.text(0.05, 0.88, 'Fold', transform=ax3.transAxes,size=12)

ax1.plot(t, p1, c=[220/255,5/255,12/255], linewidth=2.5)
ax2.plot(t, p2, c=[78/255,178/255,101/255], linewidth=2.5)
ax2.plot(t, p3, c=[220/255,5/255,12/255], linewidth=2.5)
ax3.plot(t, p4, c=[78/255,178/255,101/255], linewidth=2.5)
ax3.plot(t, p5, c=[220/255,5/255,12/255], linewidth=2.5)

ax4.set_xlim(-0.55, 0.1)
ax4.set_ylim(-0.5, 4)
ax6.set_ylim(-0.5, 4.5)
ax4.set_xlabel('External forcing', fontsize = fontsize)
ax4.set_ylabel('System state', fontsize = fontsize)
ax5.set_xlabel('External forcing', fontsize = fontsize)
ax6.set_xlabel('External forcing', fontsize = fontsize)
ax4.tick_params(axis='both', which='major', labelsize=10)
ax5.tick_params(axis='both', which='major', labelsize=10)
ax6.tick_params(axis='both', which='major', labelsize=10)
ax4.set_xticks([])
ax4.set_yticks([])
ax5.set_yticks([])
ax6.set_yticks([])

fig.tight_layout()



for i in range(n):
    ## Solve ODE with Euler method
    x[i+1] = x[i] + dt*f(x[i], p1[i], A1, s1)
    x2[i+1] = x2[i] + dt*f(x2[i], p2[i], A1, s1)
    x3[i+1] = x3[i] + dt*f(x3[i], p3[i], A1, s1)
    x4[i+1] = x4[i] + dt*f(x4[i], p4[i], A2, s2)
    x5[i+1] = x5[i] + dt*f(x5[i], p5[i], A2, s2)

xx = np.linspace(-2.5,4.5,601)

## Plotting

ax4.plot(pp, xeqplus1, 'k', linewidth=1.5)
ax4.plot(pp, xeqminus1, 'k--', linewidth=1.5)
ax4.plot(pp, np.zeros(len(pp)), 'k', linewidth=1.5)
ax4.plot(p1, x,c=[220/255,5/255,12/255],linewidth=2.5)
ax4.plot([p0, BI1],[x[0], x[0]],'k:')
ax4.plot(p1[0], x[0], c=[0.8,0.8,0.8], marker='.', MarkerSize = 16)
ax4.plot(p1[-1], x[-1], c=[220/255,5/255,12/255], marker='.', MarkerSize = 16)
ax4.fill_between([BI1,0],-0.5,4.5,color=[0.85,0.85,0.85])
ax4.text(BI1-BI1/2, 0.94, "Basin \n instability", size=labelsize, ha="center", va="center")
ax4.text(-0.4, 2.5, "Base state", size=labelsize, ha="center", va="center",rotation=15)
ax4.text(-0.38, 0.2, "Alternative state", size=labelsize, ha="center", va="center")


ax5.plot(pp, xeqplus1, 'k', linewidth=1.5)
ax5.plot(pp, xeqminus1, 'k--', linewidth=1.5)
ax5.plot(pp, np.zeros(len(pp)), 'k', linewidth=1.5)
ax5.plot(p2, x2,c=[78/255,178/255,101/255],linewidth=2.5)
ax5.plot(p3, x3,c=[220/255,5/255,12/255],linewidth=2.5)
ax5.plot([p0, BI1],[x[0], x[0]],'k:')
ax5.plot(p2[0], x2[0], c=[0.8,0.8,0.8], marker='.', MarkerSize = 16)
ax5.plot(p2[-1], x2[-1], c=[78/255,178/255,101/255], marker='.', MarkerSize = 16)
ax5.plot(p3[-1], x3[-1], c=[220/255,5/255,12/255], marker='.', MarkerSize = 16)
ax5.fill_between([BI1,0],-0.5,4.5,color=[0.85,0.85,0.85])
ax5.text(BI1+0.05, 0.94, "Basin \n instability", size=labelsize, ha="center", va="center")


ax6.plot(pp, xeqplus2, 'k', linewidth=1.5)
ax6.plot(pp, xeqminus2, 'k--', linewidth=1.5)
ax6.plot(pp, np.zeros(len(pp)), 'k', linewidth=1.5)
ax6.plot(p4, x4,c=[78/255,178/255,101/255],linewidth=2.5)
ax6.plot(p5, x5,c=[220/255,5/255,12/255],linewidth=2.5)
ax6.plot([p0+Deltap2, BI2],[x4[0], x4[0]],'k:')
ax6.plot(p4[0], x4[0], c=[0.8,0.8,0.8], marker='.', MarkerSize = 16)
ax6.plot(p4[-1], x4[-1], c=[78/255,178/255,101/255], marker='.', MarkerSize = 16)
ax6.plot(p5[-1], x5[-1], c=[220/255,5/255,12/255], marker='.', MarkerSize = 16)
ax6.fill_between([-1,BI2],-0.5,4.5,color=[0.85,0.85,0.85])
ax6.text(BI2-0.1, 1.1, "Basin \n instability", size=labelsize, ha="center", va="center")


## Adding arrows
style = "Simple, tail_width=0.5, head_width=4, head_length=8"
kw = dict(arrowstyle=style)

idx = np.nanargmin(np.abs(p1 - -0.25)+np.abs(x - 2.7))
arrow0 = p1[idx+1], x[idx+1]
arrow1 = p1[idx], x[idx]
ax4.annotate('',xytext=(arrow1),xy=(arrow0),arrowprops=dict(arrowstyle="simple", color=[220/255,5/255,12/255]),size=16)
idx = np.nanargmin(np.abs(p1 - 0.001)+np.abs(x - 1))
arrow0 = p1[idx+1], x[idx+1]
arrow1 = p1[idx], x[idx]
ax4.annotate('',xytext=(arrow1),xy=(arrow0),arrowprops=dict(arrowstyle="simple", color=[220/255,5/255,12/255]),size=16)

idx = np.nanargmin(np.abs(p2 - -0.2)+np.abs(x2 - 2.7))
arrow0 = p2[idx+1], x2[idx+1]
arrow1 = p2[idx], x2[idx]
ax5.annotate('',xytext=(arrow1),xy=(arrow0),arrowprops=dict(arrowstyle="simple", color=[78/255,178/255,101/255]),size=16)

idx = np.nanargmin(np.abs(p3 - -0.35)+np.abs(x3 - 2.25))
arrow0 = p3[idx+1], x3[idx+1]
arrow1 = p3[idx], x3[idx]
ax5.annotate('',xytext=(arrow1),xy=(arrow0),arrowprops=dict(arrowstyle="simple", color=[220/255,5/255,12/255]),size=16)
idx = np.nanargmin(np.abs(p3 - -0.25)+np.abs(x3 - 2.274))
arrow0 = p3[idx+1], x3[idx+1]
arrow1 = p3[idx], x3[idx]
ax5.annotate('',xytext=(arrow1),xy=(arrow0),arrowprops=dict(arrowstyle="simple", color=[220/255,5/255,12/255]),size=16)
idx = np.nanargmin(np.abs(p3 - -0.09)+np.abs(x3 - 1))
arrow0 = p3[idx+1], x3[idx+1]
arrow1 = p3[idx], x3[idx]
ax5.annotate('',xytext=(arrow1),xy=(arrow0),arrowprops=dict(arrowstyle="simple", color=[220/255,5/255,12/255]),size=16)

idx = np.nanargmin(np.abs(p4 - -0.4)+np.abs(x4 - 3.4))
arrow0 = p4[idx+1], x4[idx+1]
arrow1 = p4[idx], x4[idx]
ax6.annotate('',xytext=(arrow1),xy=(arrow0),arrowprops=dict(arrowstyle="simple", color=[78/255,178/255,101/255]),size=16)

idx = np.nanargmin(np.abs(p5 - -0.25)+np.abs(x5 - 2.45))
arrow0 = p5[idx+1], x5[idx+1]
arrow1 = p5[idx], x5[idx]
ax6.annotate('',xytext=(arrow1),xy=(arrow0),arrowprops=dict(arrowstyle="simple", color=[220/255,5/255,12/255]),size=16)
idx = np.nanargmin(np.abs(p5 - -0.25)+np.abs(x5 - 2.475))
arrow0 = p5[idx+1], x5[idx+1]
arrow1 = p5[idx], x5[idx]
ax6.annotate('',xytext=(arrow1),xy=(arrow0),arrowprops=dict(arrowstyle="simple", color=[220/255,5/255,12/255]),size=16)
idx = np.nanargmin(np.abs(p5 - -0.5)+np.abs(x5 - 1.5))
arrow0 = p5[idx+1], x5[idx+1]
arrow1 = p5[idx], x5[idx]
ax6.annotate('',xytext=(arrow1),xy=(arrow0),arrowprops=dict(arrowstyle="simple", color=[220/255,5/255,12/255]),size=16)

## Labelling
ax1.text(0.94, 0.93, 'a)', transform=ax1.transAxes,size=12, weight='bold')
ax2.text(0.94, 0.93, 'b)', transform=ax2.transAxes,size=12, weight='bold')
ax3.text(0.94, 0.93, 'c)', transform=ax3.transAxes,size=12, weight='bold')

ax4.text(0.94, 0.93, 'd)', transform=ax4.transAxes,size=12, weight='bold')
ax5.text(0.94, 0.93, 'e)', transform=ax5.transAxes,size=12, weight='bold')
ax6.text(0.94, 0.93, 'f)', transform=ax6.transAxes,size=12, weight='bold')


ax4.plot(0, A1, c='k', marker='.', MarkerSize = 16,zorder=1)
ax5.plot(0, A1, c='k', marker='.', MarkerSize = 16,zorder=1)
ax6.plot(0, A2, c='k', marker='.', MarkerSize = 16,zorder=1)

ax4.text(0.88, 0.8, 'Fold', transform=ax4.transAxes,size=12)
ax5.text(0.88, 0.8, 'Fold', transform=ax5.transAxes,size=12)

ax6.text(0.88, 0.4, 'Fold', transform=ax6.transAxes,size=12)