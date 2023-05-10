# -*- coding: utf-8 -*-
"""
Created on Thu May  4 09:49:46 2023

@author: Paul

Figure 5. Tipping diagrams for the Plant-Herbivore and AMOC models
"""

#import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


plt.rcParams["font.family"] = "serif"
fontsize = 12
labelsize=10


## Import P-H data
data = pd.read_excel("tipping_diagrams_AMOC_and_PH_v2.xlsx", 'PH_model')
data2=data.to_numpy()

r = np.concatenate([data2[0,1:]]).astype(None)
Delta_rho_1 = np.concatenate([data2[1,1:]]).astype(None)
Delta_rho_2 = np.concatenate([data2[2,1:]]).astype(None)

SN = data2[4,1]
BI = data2[5,1]

#### Plotting critical boundaries for P-H model
fig, ax = plt.subplots(1,2,figsize=(11,4.5))
ax[0].plot([SN,SN],[1E-5,1],'k',linewidth=0.5,label='Fold')
ax[0].plot([BI,BI],[1E-5,1],'k',linewidth=0.5,label='BI')
ax[0].plot(Delta_rho_2,r,'k')
ax[0].plot(Delta_rho_1,r,'k--')
ax[0].fill_betweenx(r,Delta_rho_1,Delta_rho_2,color=[78/255,178/255,101/255],alpha=0.4,edgecolor='none')
ax[0].fill_betweenx(r,Delta_rho_2,10,color=[220/255,5/255,12/255],alpha=0.4)
ax[0].set_xlim(0,1.5)
ax[0].set_ylim(1E-3,1E0)
ax[0].set_yscale('log')
ax[0].set_xlabel('Peak change in environmental conditions',fontsize=fontsize)
ax[0].set_ylabel('Rate parameter of environmental \n conditions [per day]',fontsize=fontsize)

ax[0].text(0.03, 0.93, 'a)', transform=ax[0].transAxes,size=12,fontsize=fontsize, weight='bold')
ax[0].text(BI-0.026, 1E-2, "Basin instability boundary", size=labelsize, rotation=90, ha="center", va="center")
ax[0].text(SN-0.026, 0.5, "Fold", size=labelsize, rotation=90, ha="center", va="center")




## Import AMOC data
data = pd.read_excel("tipping_diagrams_AMOC_and_PH_v2.xlsx", 'AMOC')
data2=data.to_numpy()

r = np.concatenate([data2[0,1:]]).astype(None)
Delta_H_1 = np.concatenate([data2[1,1:]]).astype(None)
Delta_H_2 = np.concatenate([data2[2,1:]]).astype(None)

Hopf = data2[4,1]
SN = data2[5,1]
BI = data2[6,1]

## Plotting critical boundaries for AMOC
ax[1].plot([SN,SN],[1E-6,1],'k',linewidth=0.5,label='Fold')
ax[1].plot([Hopf,Hopf],[1E-6,1],'k',linewidth=0.5,label='Hopf')
ax[1].plot([BI,BI],[1E-6,1],'k',linewidth=0.5,label='BI')
ax[1].plot(Delta_H_1,r,'k--')
ax[1].plot(Delta_H_2,r,'k')
ax[1].fill_betweenx(r,Delta_H_1,Delta_H_2,color=[78/255,178/255,101/255],alpha=0.4,edgecolor='none')
ax[1].fill_betweenx(r,Delta_H_2,10,color=[220/255,5/255,12/255],alpha=0.4)
ax[1].set_xlim(0.35,0.45)
ax[1].set_ylim(1E-6,1E-1)
ax[1].set_yscale('log')
ax[1].set_xlabel('Peak change in freshwater hosing [Sv]',fontsize=fontsize)
ax[1].set_ylabel('Rate parameter of freshwater hosing \n [per year]',fontsize=fontsize)


ax[1].text(SN-0.002, 0.03, "Fold", size=labelsize, rotation=90, ha="center", va="center")
ax[1].text(Hopf-0.002, 0.03, "Hopf", size=labelsize, rotation=90, ha="center", va="center")
ax[1].text(BI-0.002, 4.5E-5, "Basin instability boundary", size=labelsize, rotation=90, ha="center", va="center")

ax[1].text(0.03, 0.93, 'b)', transform=ax[1].transAxes,size=12,fontsize=fontsize, weight='bold')


fig.tight_layout()