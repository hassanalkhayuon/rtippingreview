# -*- coding: utf-8 -*-
"""
Created on Thu May  4 09:58:58 2023

@author: Paul

Figure 6. Tipping for ramp and return forcing profiles in the power grid model.
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


def BusModel_4D_PG_IVP(t,x,P_1,Q_shift,r,RETURN):

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
    if RETURN:
        z[4] = (-Q_shift*r*np.tanh(r*t)/np.cosh(r*t))
    else:
        z[4] = (-Q_shift*r*np.tanh(r*t)/np.cosh(r*t))*(t<0)+0*(t>=0)
    
    return(z)

def Forcing(t, Q1, Q_shift, r):
    """
    External forcing ramp
    """    
    Q = (Q1 + Q_shift/np.cosh(r*t))*(t<=0) + (Q1+Q_shift)*(t>0) 
    return (Q)

def Forcing2(t, Q1, Q_shift, r):
    """
    External forcing return
    """    
    Q = (Q1 + Q_shift/np.cosh(r*t))
    return (Q)

def event(t,x,P_1,Q_shift,r):
    return x[3]
    

## System parameter
P_1 = 5
## Initial variable values
delta_m = 0.1119
w_m = 0
delta = -0.097
V = 1.5149
Q_1 = 5

## Forcing parameters
Q_shift = 5
r1 = 50
r2 = 130
r3 = 400


## Set up for initial value problem solver
tspan = [-1,1]
h = 0.0001
t = np.arange(tspan[0],tspan[1]+1,h)
tol = 1E-4

colours2 = ['b',(0.6,0,0.6),'r']


sol1 = solve_ivp(BusModel_4D_PG_IVP,tspan,[delta_m, w_m, delta, V, Q_1],args=(P_1,Q_shift,r1,False),teval=t,max_step=h)
sol2 = solve_ivp(BusModel_4D_PG_IVP,tspan,[delta_m, w_m, delta, V, Q_1],args=(P_1,Q_shift,r1,True),teval=t,max_step=h)
sol3 = solve_ivp(BusModel_4D_PG_IVP,tspan,[delta_m, w_m, delta, V, Q_1],args=(P_1,Q_shift,r2,False),teval=t,max_step=h)
sol4 = solve_ivp(BusModel_4D_PG_IVP,tspan,[delta_m, w_m, delta, V, Q_1],args=(P_1,Q_shift,r2,True),teval=t,max_step=h)
sol5 = solve_ivp(BusModel_4D_PG_IVP,tspan,[delta_m, w_m, delta, V, Q_1],args=(P_1,Q_shift,r3,False),teval=t,max_step=h)
sol6 = solve_ivp(BusModel_4D_PG_IVP,tspan,[delta_m, w_m, delta, V, Q_1],args=(P_1,Q_shift,r3,True),teval=t,max_step=h)
    
## Initialise figure
fig=plt.figure(figsize=(11,4.5))
gs=GridSpec(3,2,hspace=0,height_ratios=[1,0.1,1])

ax1=fig.add_subplot(gs[0,0]) 
ax2=fig.add_subplot(gs[2,0])
ax3=fig.add_subplot(gs[:,1])

## Plotting
ax1.plot(t[t<0]+0.15,Forcing(t[t<0], Q_1, Q_shift, r1),c=[25/255,101/255,176/255],ls='--')
ax1.plot(t[t<0]+0.15,Forcing(t[t<0], Q_1, Q_shift, r2),c=[238/255,128/255,38/255],ls='--')
ax1.plot(t[t<0]+0.15,Forcing(t[t<0], Q_1, Q_shift, r3),c=[174/255,118/255,163/255],ls='--')
ax1.plot(t[t>=0]+0.15,Forcing(t[t>=0], Q_1, Q_shift, r1),c='k',ls='--')
ax1.plot(t+0.15,Forcing2(t, Q_1, Q_shift, r1),c=[25/255,101/255,176/255])
ax1.plot(t+0.15,Forcing2(t, Q_1, Q_shift, r2),c=[238/255,128/255,38/255])
ax1.plot(t+0.15,Forcing2(t, Q_1, Q_shift, r3),c=[174/255,118/255,163/255])
ax1.set_ylabel('Power demand',fontsize=fontsize,labelpad=10)
ax1.set_xticklabels([])
ax1.set_xlim(0,0.55)

ax2.plot(sol1.t+0.15,sol1.y[3,:],c=[25/255,101/255,176/255],ls='--')
ax2.plot(sol3.t+0.15,sol3.y[3,:],c=[238/255,128/255,38/255],ls='--')
ax2.plot(sol5.t+0.15,sol5.y[3,:],c=[174/255,118/255,163/255],ls='--')
ax2.plot(sol2.t+0.15,sol2.y[3,:],c=[25/255,101/255,176/255])
ax2.plot(sol4.t+0.15,sol4.y[3,:],c=[238/255,128/255,38/255])
ax2.plot(sol6.t+0.15,sol6.y[3,:],c=[174/255,118/255,163/255])
ax2.set_xlabel('Time',fontsize=fontsize,labelpad=10)
ax2.set_ylabel('Voltage',fontsize=fontsize,labelpad=10)
ax2.set_xlim(0,0.55) 
ax2.set_ylim(0,2) 



## Import critical boundary data
mat1 = loadmat('Power_grid_model_monotone_forcing_tracking_P5_data.mat')
mat2 = loadmat('Power_grid_model_nonmonotone_forcing_tracking_P5_data.mat')
mat3 = loadmat('Power_grid_model_monotone_forcing_tipping_P5_data.mat')
mat4 = loadmat('Power_grid_model_nonmonotone_forcing_tipping_P5_data.mat')

mat5 = loadmat('Power_grid_model_monotone_forcing_tipping_slow_rates_P5_data.mat')
mat6 = loadmat('Power_grid_model_nonmonotone_forcing_tipping_slow_rates_P5_data.mat')
mat7 = loadmat('Power_grid_model_monotone_forcing_tracking_slow_rates_P5_data.mat')
mat8 = loadmat('Power_grid_model_nonmonotone_forcing_tracking_slow_rates_P5_data.mat')


Q_shifts_tip = np.concatenate((mat5['Q_shifts'][0,:],mat3['Q_shifts'][0,34:]))
rs_tip = np.concatenate((mat5['rs'][0,:],mat3['rs'][0,34:]))
Q_returns_tip = np.concatenate((mat6['Q_shifts'][0,:],mat4['Q_shifts'][0,34:]))
rs_returns_tip = np.concatenate((mat6['rs'][0,:],mat4['rs'][0,34:]))
Q_shifts_track = np.concatenate((mat7['Q_shifts'][0,:],mat1['Q_shifts'][0,34:]))
rs_track = np.concatenate((mat7['rs'][0,:],mat1['rs'][0,34:]))
Q_returns_track = np.concatenate((mat8['Q_shifts'][0,:],mat2['Q_shifts'][0,34:]))
rs_returns_track = np.concatenate((mat8['rs'][0,:],mat2['rs'][0,34:]))

f = interpolate.interp1d(rs_track, Q_shifts_track)
f2 = interpolate.interp1d(rs_returns_track, Q_returns_track)
f3 = interpolate.interp1d(rs_tip, Q_shifts_tip)
f4 = interpolate.interp1d(rs_returns_tip, Q_returns_tip)

x = np.geomspace(1,1000,10001)
y1 = f(x)
y2 = f2(x)
y3 = f3(x)
y4 = f4(x)

mat9 = loadmat('Power_grid_model_nonmonotone_forcing_tipping_nonsmooth_int_lower_curve_P5_data.mat')
mat10 = loadmat('Power_grid_model_nonmonotone_forcing_tipping_nonsmooth_int_upper_curve_P5_data.mat')


Q_shifts_nonmon_tip_lower = mat9['Q_shifts'][0,:]
rs_nonmon_tip_lower = mat9['rs'][0,:]
Q_shifts_nonmon_tip_upper = mat10['Q_shifts'][0,:]
rs_nonmon_tip_upper = mat10['rs'][0,:]

f5 = interpolate.interp1d(Q_shifts_nonmon_tip_lower, rs_nonmon_tip_lower)
f6 = interpolate.interp1d(Q_shifts_nonmon_tip_upper, rs_nonmon_tip_upper)

x2 = np.linspace(4.2701,6.29,10001)
y5 = f5(x2)
y6 = f6(x2)

idx=~np.isnan(y5)

X = np.concatenate((y4[x<50],x2[idx][::-1],x2,y4[x>100]))
Y = np.concatenate((x[x<50],y5[idx][::-1],y6,x[x>100]))

## Plotting critical boundaries

ax3.plot(y1,x,c='k',linestyle='dashed')
ax3.plot(y2,x,c='k')
ax3.plot(X,Y,c='k')
ax3.plot([5.34,5.34],[1,1000],'k',lw=0.5)
BI = 3.52880859
ax3.plot([BI, BI],[1,1000],'k',lw=0.5)

ax3.fill_betweenx(x,y2,y1,color=[78/255,178/255,101/255],edgecolor='none',alpha=0.4)
ax3.fill_betweenx(x[x<50],y2[x<50],y4[x<50],color=[220/255,5/255,12/255],edgecolor='none',alpha=0.4)
ax3.fill_betweenx(x[x>100],y2[x>100],y4[x>100],color=[220/255,5/255,12/255],edgecolor='none',alpha=0.4)
ax3.fill_between(x2,y6,np.nanmin(x[x>100]),color=[220/255,5/255,12/255],edgecolor='none',alpha=0.4)
ax3.fill_betweenx(Y,X,7.5,color=[220/255,5/255,12/255],edgecolor='none',alpha=0.7)
ax3.set_xlim(3,7)
ax3.set_ylim(1,1000)
ax3.set_yscale('log')
ax3.set_xlabel('Peak change in power demand',fontsize=fontsize,labelpad=10)
ax3.set_ylabel('Rate parameter of power demand',fontsize=fontsize,labelpad=10)


sns.despine()
ax3.spines["top"].set_visible(True)
ax3.spines["right"].set_visible(True)

ax3.plot(Q_shift,r1,Marker='o',ms=10,c=[25/255,101/255,176/255])
ax3.plot(Q_shift,r2,Marker='o',ms=10,c=[238/255,128/255,38/255])
ax3.plot(Q_shift,r3,Marker='o',ms=10,c=[174/255,118/255,163/255])
ax3.text(5.25, 600, "Fold", size=labelsize, rotation=90, ha="center", va="center")
ax3.text(3.44, 10, "Basin instability boundary", size=labelsize, rotation=90, ha="center", va="center")

ax1.text(0.03, 0.86, 'a)', transform=ax1.transAxes,size=12,fontsize=fontsize, weight='bold')
ax2.text(0.03, 0.86, 'b)', transform=ax2.transAxes,size=12,fontsize=fontsize, weight='bold')
ax3.text(0.03, 0.93, 'c)', transform=ax3.transAxes,size=12,fontsize=fontsize, weight='bold')

fig.tight_layout()