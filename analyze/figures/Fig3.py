#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 1 14:45:51 2021

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


#Load stimulation Data
file_shapes=np.load('../collectedData/CollectResults_shapesLong.npz')
collect_Iso_shapes=file_shapes['collect_Iso']
collect_Ispindles_shapes=file_shapes['collect_Ispindles']
collectr_Iso_shapes=file_shapes['collectr_Iso']
collectr_Ispindles_shapes=file_shapes['collectr_Ispindles']
collect_pc_shapes=file_shapes['collect_pc']
collect_pcp_shapes=file_shapes['collect_pcp']
collect_pso_shapes=file_shapes['collect_pso']
collect_psp_shapes=file_shapes['collect_psp']
collect_coincidences_shapes=file_shapes['collect_coincidences']
collect_so_shapes=file_shapes['collect_so']
collect_spindles_shapes=file_shapes['collect_spindles']
collectr_pc_shapes=file_shapes['collectr_pc']
collectr_pcp_shapes=file_shapes['collectr_pcp']
collectr_pso_shapes=file_shapes['collectr_pso']
collectr_psp_shapes=file_shapes['collectr_psp']
collectr_coincidences_shapes=file_shapes['collectr_coincidences']
collectr_so_shapes=file_shapes['collectr_so']
collectr_spindles_shapes=file_shapes['collectr_spindles']
y_axis_shapes=[0,1,2,3,4,5]
labels_shapes=['DR','HT','Gauss','Rect','Tri','RR']


#Figure Layout
fig=plt.figure(figsize=(5.2,3.7))
gs=gridspec.GridSpec(3,6,figure=fig,width_ratios=[1,1,1,1,1,1],height_ratios=[1,0.15,1],wspace=1.15, hspace=0.4)
ax11=fig.add_subplot(gs[1,:])
axA1=fig.add_subplot(gs[0,0:3])
axA2=fig.add_subplot(gs[0,3:6])
axB1=fig.add_subplot(gs[2,0:2])
axC1=fig.add_subplot(gs[2,2:4])
axD1=fig.add_subplot(gs[2,4:6])


#Calculate shapes
A=np.array([10,17.32,17.32,17.32,15.897,15.885,14.883,14.142,16.330,14.142,20.186,30.001,12.248,12.248])
A=A*2
duration=0.1
fs=1000

sigma=duration**2/10 #5 sigma width
#time array
t=np.linspace(0,duration,int(duration*fs+1))
t1=np.linspace(0,duration+2/fs,int(duration*fs+3))
zeroshape=np.zeros((np.shape(t)[0]+2,))
first_half=np.zeros_like(t)
second_half=np.zeros_like(t)
first_quarter=np.zeros_like(t)
last_quarter=np.zeros_like(t)
middle_half=np.zeros_like(t)
first_half[0:int(duration*fs//2)+1]=1
second_half[int(duration*fs//2+1):int(duration*fs+1)]=t[int(duration*fs//2)+1:int(duration*fs+1)]-duration/2
first_quarter[0:int(duration*fs//4+1)]=t[0:int(duration*fs//4+1)]
middle_half[int(duration*fs//4+1):int(duration*fs+1)]=1
last_quarter[int(duration*fs*3//4+1):int(duration*fs+1)]=t[int(duration*fs*3//4+1):int(duration*fs+1)]-3*duration/4
#shapes
rectangular=np.ones_like(t)*A[0]
triangular=A[1]-2/duration*A[1]*np.abs(t-duration/2)
increase_ramp=A[2]*t/duration
decrease_ramp=A[3]*(1-t/duration)
half_trapece=np.ones_like(t)*A[12]-A[12]*(2*second_half/duration)
gaussian=A[4]*np.exp(-(t-duration/2)**2/sigma)
# Plot Shapes
zeroshape[1:-1]=decrease_ramp
axA1.plot(t1,zeroshape)
zeroshape[1:-1]=half_trapece
axA1.plot(t1,zeroshape)
zeroshape[1:-1]=gaussian
axA1.plot(t1,zeroshape)
zeroshape[1:-1]=rectangular
axA2.plot(t1,zeroshape,color=cm.tab10(3))
zeroshape[1:-1]=triangular
axA2.plot(t1,zeroshape,color=cm.tab10(4))
zeroshape[1:-1]=increase_ramp
axA2.plot(t1,zeroshape,color=cm.tab10(5))
axA1.set_ylim([0,40])
axA1.set_xlim([-0.01,0.12])
axA1.set_ylabel('Pulse Amplitude $s^{-1}$',fontsize=8,labelpad=0)
axA1.tick_params(axis='both', which='major', labelsize=8)
axA1.set_xlabel('Time (s)',fontsize=8,labelpad=-0.1)
axA1.text(-0.24,0.98,'A',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axA1.transAxes)

axA2.set_ylim([0,40])
axA2.set_xlim([-0.01,0.12])
# axA2.set_ylabel('Pulse Amplitude $s^{-1}$',fontsize=8,labelpad=0)
axA2.tick_params(axis='both', which='major', labelsize=8)
axA2.set_xlabel('Time (s)',fontsize=8,labelpad=-0.1)
axA2.text(-0.18,0.98,'B',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axA2.transAxes)


axB1=sp.bivariable_std(collect_Iso_shapes,collect_so_shapes,axB1,stdx_zero=stdSOb,meany_zero=np.mean(n_SO_white),stdy_zero=np.std(n_SO_white),marker='s',cmap=plt.cm.tab10 ,flagstdx=True,flagstdy=True,relative_diference=False,labels=['1','2','3','4','5','6'],meanzero_type='line_inf',colormap='direct')
axB1.set_ylabel('$N_{SO}$/min.',fontsize=8,labelpad=0)
axB1.set_xlabel('$I^{(SO)}$ (a. u.)',fontsize=8,labelpad=0)
axB1.set_ylim([10,25])
axB1.set_yticks([10,15,20,25])
axB1.set_xlim([-0.05,0.2])
axB1.text(-0.35,1,'C',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axB1.transAxes)
axB1.tick_params(axis='both', which='major', labelsize=8)
hB1,lB1=axB1.get_legend_handles_labels()
hB,lB=axB1.get_legend_handles_labels()

axC1=sp.bivariable_std(collect_Ispindles_shapes,collect_spindles_shapes,axC1,stdx_zero=stdSpindlesb,meany_zero=np.mean(n_spindles_white),stdy_zero=np.std(n_spindles_white),marker='s',cmap=plt.cm.tab10 ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_shapes,meanzero_type='line_inf',colormap='direct')
axC1.set_ylabel('$N_{SP}$/min.',fontsize=8,labelpad=0)
axC1.set_xlabel('$I^{(SP)}$ (a. u.)',fontsize=8,labelpad=0)
axC1.set_ylim([4,10])
axC1.set_xlim([-0.01,0.1])
axC1.set_yticks([4,6,8,10])
axC1.tick_params(axis='both', which='major', labelsize=8)
axC1.text(-0.3,1,'D',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axC1.transAxes)

axD1=sp.bivariable_std(collect_pso_shapes,collect_pcp_shapes,axD1,meanx_zero=np.mean(P_SO_white),stdx_zero=np.std(P_SO_white),meany_zero=np.mean(P_PCP_white),stdy_zero=np.std(P_PCP_white),marker='s',cmap=plt.cm.tab10 ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_shapes,meanzero_type='line_inf',colormap='direct')
axD1.set_ylabel('$P(C|SP)$',fontsize=8,labelpad=0)
axD1.set_xlabel('$P(SO)$',fontsize=8,labelpad=0)
axD1.set_ylim([0.1,0.5])
axD1.set_xlim([0.1,0.4])
axD1.plot([0,1],[0,1],color='gray')
axD1.set_yticks([0.1,0.2,0.3,0.4,0.5])
axD1.set_xticks([0.1,0.2,0.3,0.4])
axD1.tick_params(axis='both', which='major', labelsize=8)
axD1.text(-0.32,1,'E',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axD1.transAxes)

ax11.set_axis_off()
legfig1=ax11.legend(hB1,lB1,fontsize=8,ncol=7,loc='upper left',bbox_to_anchor=(-0.01,0.02,1,1),columnspacing=1.67)


fig.savefig('./output/Fig3.eps',dpi=300,bbox_inches='tight',bbox_extra_artists=[legfig1])
