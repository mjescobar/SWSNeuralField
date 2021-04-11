#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 14:52:18 2021

@author: felipe
"""
import sys, os
sys.path.append(os.path.abspath(os.path.join('..', 'plotting')))
from cycler import cycler
import numpy as np
import scipy.stats as stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import statistic_plots as sp
import regression as regression
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
rc('text', usetex=True)
timeSim=15.166
data_length=90999

#Load Baseline Data
filewhite=np.load('../collectedData/SWSBaselinesLong.npz') 
n_spindles_whiteall=filewhite['n_spindles_all'][0]/timeSim
n_Delta_whiteall=filewhite['n_Delta_all'][0]/timeSim
n_SO_whiteall=filewhite['n_SO_all'][0]/timeSim
n_spindles_white=filewhite['n_spindles']/timeSim
n_Delta_white=filewhite['n_Delta']/timeSim
n_SO_white=filewhite['n_SO']/timeSim
n_spindlesn_white=filewhite['n_spindlesn'][0]/timeSim
n_Deltan_white=filewhite['n_Deltan']/timeSim
coincident_centers_white=filewhite['coincident_centers']
SO_Dpeak_phase_white=filewhite['SO_Dpeak_phase'][0,:]
SpindlesSOcount_white=filewhite['SpindlesSOcount']
SpindlesSOcount_whiteall=filewhite['SpindlesSOcount_all'][0]
Spindles_centers_phase_white=filewhite['Spindles_centers_phase'][0,:]
Spindles_coincenters_phase_white=filewhite['Spindles_coincenters_phase'][0,:]
Spindles_coincenters_time_white=filewhite['Spindles_coincenters_time'][0,:]
Spindles_coincidences_phase_white=filewhite['Spindles_coincidences_phase'][0,:]
coincidences_white=filewhite['coincidences']
Spindles_centers_white=filewhite['Spindles_centers'][0,:]
P_C_white=filewhite['P_C'][0,:]
P_SP_white=filewhite['P_SP'][0,:]
P_SO_white=filewhite['P_SO'][0,:]
P_PCP_white=SpindlesSOcount_white/(n_spindles_white*timeSim)
fileSO=np.load('../collectedData/SWSbaselinesSOlong.npz')
SO_centers_white=fileSO['SO_centers_white']
SOarray_white=fileSO['SOarray_white']
c_white=int(fileSO['c_white'])
m_white=int(fileSO['m_white'])
tSO=np.arange(-2,2,0.01)
#Calculate Baseline statistics
baselineData=np.load('../collectedData/baselineDatalong.npz')
sWavelet=baselineData['sWavelet']
freqsWavelet=baselineData['freqsWavelet']
freqdot5=np.argwhere(freqsWavelet>=0.5)[-1][0]
freqone=np.argwhere(freqsWavelet>=1.25)[-1][0]
freqeleven=np.argwhere(freqsWavelet>=9)[-1][0]
freqsixteen=np.argwhere(freqsWavelet>=16)[-1][0]
stdSOb=np.std(np.sum(sWavelet[:,freqone:freqdot5],axis=1)/np.mean(np.sum(sWavelet[:,freqone:freqdot5],axis=1)))
stdSpindlesb=np.std(np.sum(sWavelet[:,freqsixteen:freqeleven+1],axis=1)/np.mean(np.sum(sWavelet[:,freqsixteen:freqeleven+1],axis=1)))

#Load Data
file_energy=np.load('../collectedData/CollectResults_energyLong.npz')
collect_Iso_energy=file_energy['collect_Iso']
collect_Ispindles_energy=file_energy['collect_Ispindles']
collectr_Iso_energy=file_energy['collectr_Iso']
collectr_Ispindles_energy=file_energy['collectr_Ispindles']
collect_pc_energy=file_energy['collect_pc']
collect_pcp_energy=file_energy['collect_pcp']
collect_pso_energy=file_energy['collect_pso']
collect_psp_energy=file_energy['collect_psp']
collect_coincidences_energy=file_energy['collect_coincidences']
collect_so_energy=file_energy['collect_so']
collect_spindles_energy=file_energy['collect_spindles']
collectr_pc_energy=file_energy['collectr_pc']
collectr_pcp_energy=file_energy['collectr_pcp']
collectr_pso_energy=file_energy['collectr_pso']
collectr_psp_energy=file_energy['collectr_psp']
collectr_coincidences_energy=file_energy['collectr_coincidences']
collectr_so_energy=file_energy['collectr_so']
collectr_spindles_energy=file_energy['collectr_spindles']
y_axis_energy=np.array([10,20,40,60,80,100])  
labels_energy=['10','20','40','60','80','100']  


file_duration=np.load('../collectedData/CollectResults_durationLong.npz')
collect_Iso_duration=file_duration['collect_Iso']
collect_Ispindles_duration=file_duration['collect_Ispindles']
collectr_Iso_duration=file_duration['collectr_Iso']
collectr_Ispindles_duration=file_duration['collectr_Ispindles']
collect_pc_duration=file_duration['collect_pc']
collect_pcp_duration=file_duration['collect_pcp']
collect_pso_duration=file_duration['collect_pso']
collect_psp_duration=file_duration['collect_psp']
collect_coincidences_duration=file_duration['collect_coincidences']
collect_so_duration=file_duration['collect_so']
collect_spindles_duration=file_duration['collect_spindles']
collectr_pc_duration=file_duration['collectr_pc']
collectr_pcp_duration=file_duration['collectr_pcp']
collectr_pso_duration=file_duration['collectr_pso']
collectr_psp_duration=file_duration['collectr_psp']
collectr_coincidences_duration=file_duration['collectr_coincidences']
collectr_so_duration=file_duration['collectr_so']
collectr_spindles_duration=file_duration['collectr_spindles']
y_axis_duration=np.array([0.05,0.075,0.15,0.30])
labels_duration=['0.05 s','0.07 s','0.15 s','0.30 s']


#Figure Layout
fig1=plt.figure(figsize=(5.2,4))
gs1=gridspec.GridSpec(4,3,figure=fig1,width_ratios=[1,1,1],height_ratios=[1,0.14,1,0.14])
gs1.update(wspace=0.4, hspace=0.4) # set the spacing between axes.

ax33=fig1.add_subplot(gs1[1,:])
axB3=fig1.add_subplot(gs1[0,0])
axC3=fig1.add_subplot(gs1[0,1])
axD3=fig1.add_subplot(gs1[0,2])
ax22=fig1.add_subplot(gs1[3,:])
axB2=fig1.add_subplot(gs1[2,0])
axC2=fig1.add_subplot(gs1[2,1])
axD2=fig1.add_subplot(gs1[2,2])

axB2=sp.bivariable_std(collect_Iso_energy,collect_so_energy,axB2,stdx_zero=stdSOb,meany_zero=np.mean(n_SO_white),stdy_zero=np.std(n_SO_white),marker='s',cmap=plt.cm.Blues ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_energy,meanzero_type='line_inf')
axB2.set_ylabel('$N_{SO}$/min.',fontsize=8,labelpad=0)
axB2.set_xlabel('$I^{(SO)}$ (a. u.)',fontsize=8,labelpad=0)
axB2.set_ylim([10,25])
axB2.set_xlim([-0.05,0.6])
axB2.tick_params(axis='both', which='major', labelsize=8)
axB2.text(-0.33,0.93,'D',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axB2.transAxes)
hB,lB=axB2.get_legend_handles_labels()
hB2,lB2=axB2.get_legend_handles_labels()
hB2.remove(hB2[0])
lB2.remove(lB2[0])
for h in hB2:
    hB.append(h)
for l in lB2:
    lB.append(l)

axC2=sp.bivariable_std(collect_Ispindles_energy,collect_spindles_energy,axC2,stdx_zero=stdSpindlesb,meany_zero=np.mean(n_spindles_white),stdy_zero=np.std(n_spindles_white),marker='s',cmap=plt.cm.Blues ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_energy,meanzero_type='line_inf')
axC2.set_ylabel('$N_{SP}$/min.',fontsize=8,labelpad=0)
axC2.set_xlabel('$I^{(SP)}$ (a. u.)',fontsize=8,labelpad=0)
axC2.set_ylim([5,15])
axC2.set_xlim([-0.01,0.3])
axC2.set_yticks([5,10,15])
axC2.tick_params(axis='both', which='major', labelsize=8)
axC2.text(-0.33,0.93,'E',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axC2.transAxes)


axD2=sp.bivariable_std(collect_pso_energy,collect_pcp_energy,axD2,meanx_zero=np.mean(P_SO_white),stdx_zero=np.std(P_SO_white),meany_zero=np.mean(P_PCP_white),stdy_zero=np.std(P_PCP_white),marker='s',cmap=plt.cm.Blues ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_energy,meanzero_type='line_inf')
axD2.set_ylabel('$P(C|SP)$',fontsize=8,labelpad=0)
axD2.set_xlabel('P(SO)',fontsize=8,labelpad=0.5)
axD2.set_ylim([0.1,0.6])
axD2.set_xlim([0.1,0.4])
axD2.plot([0,1],[0,1],color='gray')
axD2.set_yticks([0.1,0.2,0.3,0.4,0.5,0.6])
axD2.set_xticks([0.1,0.2,0.3,0.4])
axD2.tick_params(axis='both', which='major', labelsize=8)
axD2.text(-0.34,0.93,'F',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axD2.transAxes)


axB3=sp.bivariable_std(collect_Iso_duration,collect_so_duration,axB3,stdx_zero=stdSOb,meany_zero=np.mean(n_SO_white),stdy_zero=np.std(n_SO_white),marker='s',cmap=plt.cm.Reds ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_duration,meanzero_type='line_inf')
axB3.set_ylabel('$N_{SO}$/min.',fontsize=8,labelpad=0)

axB3.set_ylim([10,25])
axB3.set_xlim([-0.05,0.2])
axB3.tick_params(axis='both', which='major', labelsize=8)
axB3.text(-0.33,0.93,'A',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axB3.transAxes)
hB3,lB3=axB3.get_legend_handles_labels()
hB3.remove(hB3[0])
lB3.remove(lB3[0])
for h in hB3:
    hB.append(h)
for l in lB3:
    lB.append(l)

axC3=sp.bivariable_std(collect_Ispindles_duration,collect_spindles_duration,axC3,stdx_zero=stdSpindlesb,meany_zero=np.mean(n_spindles_white),stdy_zero=np.std(n_spindles_white),marker='s',cmap=plt.cm.Reds ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_duration,meanzero_type='line_inf')
axC3.set_ylabel('$N_{SP}$/min.',fontsize=8,labelpad=0)
axC3.set_ylim([4,10])
axC3.set_xlim([-0.01,0.1])
axC3.set_yticks([4,6,8,10])
axC3.tick_params(axis='both', which='major', labelsize=8)
axC3.text(-0.33,0.93,'B',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axC3.transAxes)

axD3=sp.bivariable_std(collect_pso_duration,collect_pcp_duration,axD3,meanx_zero=np.mean(P_SO_white),stdx_zero=np.std(P_SO_white),meany_zero=np.mean(P_PCP_white),stdy_zero=np.std(P_PCP_white),marker='s',cmap=plt.cm.Reds ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_duration,meanzero_type='line_inf')
axD3.set_ylabel('$P(C|SP)$',fontsize=8,labelpad=0)
axD3.set_ylim([0.1,0.5])
axD3.set_xlim([0.1,0.4])
axD3.plot([0,1],[0,1],color='gray')
axD3.set_yticks([0.1,0.2,0.3,0.4,0.5])
axD3.set_xticks([0.1,0.2,0.3,0.4])
axD3.tick_params(axis='both', which='major', labelsize=8)
axD3.text(-0.34,0.93,'C',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axD3.transAxes)

ax22.set_axis_off()
ax22.legend(hB2,lB2,fontsize=8,ncol=7,loc='upper left',bbox_to_anchor=(-0.01,-0.7,1,1),columnspacing=2.192)

ax33.set_axis_off()
ax33.legend(hB3,lB3,fontsize=8,ncol=7,loc='upper left',bbox_to_anchor=(-0.01,0.02,1,1),columnspacing=4.382)

#Save figure
fig1.savefig('./output/Fig4.eps',dpi=300,bbox_inches='tight')