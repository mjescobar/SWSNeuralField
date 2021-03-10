#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  6 02:09:54 2021

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

#Load Data
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

baselineData=np.load('../collectedData/baselineDatalong.npz')
sWavelet=baselineData['sWavelet']
freqsWavelet=baselineData['freqsWavelet']
freqdot5=np.argwhere(freqsWavelet>=0.5)[-1][0]
freqone=np.argwhere(freqsWavelet>=1.25)[-1][0]
freqeleven=np.argwhere(freqsWavelet>=9)[-1][0]
freqsixteen=np.argwhere(freqsWavelet>=16)[-1][0]
stdSOb=np.std(np.sum(sWavelet[:,freqone:freqdot5],axis=1)/np.mean(np.sum(sWavelet[:,freqone:freqdot5],axis=1)))
stdSpindlesb=np.std(np.sum(sWavelet[:,freqsixteen:freqeleven+1],axis=1)/np.mean(np.sum(sWavelet[:,freqsixteen:freqeleven+1],axis=1)))


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
fig2=plt.figure(figsize=(5.2,6))
gs2=gridspec.GridSpec(6,3,figure=fig2,width_ratios=[1,1,1],height_ratios=[1,0.14,1,0.14,1,0.14])
gs2.update(wspace=0.5, hspace=0.4) # set the spacing between axes.
ax11r=fig2.add_subplot(gs2[1,:])
axB4=fig2.add_subplot(gs2[0,0])
axC4=fig2.add_subplot(gs2[0,1])
axD4=fig2.add_subplot(gs2[0,2])
ax33r=fig2.add_subplot(gs2[3,:])
axB6=fig2.add_subplot(gs2[2,0])
axC6=fig2.add_subplot(gs2[2,1])
axD6=fig2.add_subplot(gs2[2,2])
ax22r=fig2.add_subplot(gs2[5,:])
axB5=fig2.add_subplot(gs2[4,0])
axC5=fig2.add_subplot(gs2[4,1])
axD5=fig2.add_subplot(gs2[4,2])


axB4=sp.bivariable_std(collectr_Iso_shapes,collectr_so_shapes,axB4,stdx_zero=stdSOb,meany_zero=np.mean(n_SO_white),stdy_zero=np.std(n_SO_white),marker='s',cmap=plt.cm.tab10 ,flagstdx=True,flagstdy=True,relative_diference=False,labels=['1','2','3','4','5','6'],meanzero_type='line_inf',colormap='direct')
axB4.set_ylabel('$N_{SO}$/min.',fontsize=8,labelpad=0)
axB4.set_ylim([10,25])
axB4.set_xlim([-0.05,0.2])
axB4.tick_params(axis='both', which='major', labelsize=8)
axB4.text(-0.35,0.9,'A',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axB4.transAxes)
hB1,lB1=axB4.get_legend_handles_labels()
hB1.remove(hB1[0])
lB1.remove(lB1[0])

axC4=sp.bivariable_std(collectr_Ispindles_shapes,collectr_spindles_shapes,axC4,stdx_zero=stdSpindlesb,meany_zero=np.mean(n_spindles_white),stdy_zero=np.std(n_spindles_white),marker='s',cmap=plt.cm.tab10 ,flagstdx=True,flagstdy=True,relative_diference=False,labels=['1','2','3','4','5','6'],meanzero_type='line_inf',colormap='direct')
axC4.set_ylabel('$N_{SP}$/min.',fontsize=8,labelpad=0)
axC4.set_ylim([4,10])
axC4.set_xlim([-0.01,0.1])
axC4.set_yticks([4,6,8,10])
axC4.tick_params(axis='both', which='major', labelsize=8)
axC4.text(-0.4,0.93,'B',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axC4.transAxes)

axD4=sp.bivariable_std(collectr_pso_shapes,collectr_pcp_shapes,axD4,meanx_zero=np.mean(P_SO_white),stdx_zero=np.std(P_SO_white),meany_zero=np.mean(P_PCP_white),stdy_zero=np.std(P_PCP_white),marker='s',cmap=plt.cm.tab10 ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_shapes,meanzero_type='line_inf',colormap='direct')
axD4.set_ylabel('$P(C|SP)$',fontsize=8,labelpad=0)
axD4.set_ylim([0.1,0.5])
axD4.set_xlim([0.1,0.4])
axD4.plot([0,1],[0,1],color='gray')
axD4.set_xticks([0.1,0.2,0.3,0.4])
axD4.set_yticks([0.1,0.2,0.3,0.4,0.5])
axD4.tick_params(axis='both', which='major', labelsize=8)
axD4.text(-0.36,0.93,'C',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axD4.transAxes)

axB5=sp.bivariable_std(collectr_Iso_energy,collectr_so_energy,axB5,stdx_zero=stdSOb,meany_zero=np.mean(n_SO_white),stdy_zero=np.std(n_SO_white),marker='s',cmap=plt.cm.Blues ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_energy,meanzero_type='line_inf')
axB5.set_ylabel('$N_{SO}$/min.',fontsize=8,labelpad=0)
axB5.set_xlabel('$I^{(SO)}$ (a. u.)',fontsize=8,labelpad=0)
axB5.set_ylim([10,25])
axB5.set_xlim([-0.05,0.6])
axB5.tick_params(axis='both', which='major', labelsize=8)
axB5.text(-0.35,0.93,'G',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axB5.transAxes)
hB2,lB2=axB5.get_legend_handles_labels()
hB2.remove(hB2[0])
lB2.remove(lB2[0])

axC5=sp.bivariable_std(collectr_Ispindles_energy,collectr_spindles_energy,axC5,stdx_zero=stdSpindlesb,meany_zero=np.mean(n_spindles_white),stdy_zero=np.std(n_spindles_white),marker='s',cmap=plt.cm.Blues ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_energy,meanzero_type='line_inf')
axC5.set_ylabel('$N_{SP}$/min.',fontsize=8,labelpad=0)
axC5.set_xlabel('$I^{(SP)}$ (a. u.)',fontsize=8,labelpad=0)
axC5.set_ylim([5,15])
axC5.set_xlim([-0.01,0.3])
axC5.tick_params(axis='both', which='major', labelsize=8)
axC5.text(-0.4,0.93,'H',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axC5.transAxes)

axD5=sp.bivariable_std(collectr_pso_energy,collectr_pcp_energy,axD5,meanx_zero=np.mean(P_SO_white),stdx_zero=np.std(P_SO_white),meany_zero=np.mean(P_PCP_white),stdy_zero=np.std(P_PCP_white),marker='s',cmap=plt.cm.Blues ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_energy,meanzero_type='line_inf')
axD5.set_ylabel('$P(C|SP)$',fontsize=8,labelpad=0)
axD5.set_xlabel('P(SO)',fontsize=8,labelpad=0)
axD5.set_ylim([0.1,0.6])
axD5.set_xlim([0.1,0.4])
axD5.plot([0,1],[0,1],color='gray')
axD5.set_xticks([0.1,0.2,0.3,0.4])
axD5.set_yticks([0.1,0.2,0.3,0.4,0.5,0.6])
axD5.tick_params(axis='both', which='major', labelsize=8)
axD5.text(-0.36,0.93,'I',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axD5.transAxes)


axB6=sp.bivariable_std(collectr_Iso_duration,collectr_so_duration,axB6,stdx_zero=stdSOb,meany_zero=np.mean(n_SO_white),stdy_zero=np.std(n_SO_white),marker='s',cmap=plt.cm.Reds ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_duration,meanzero_type='line_inf')
axB6.set_ylabel('$N_{SO}$/min.',fontsize=8,labelpad=0)
axB6.set_ylim([10,25])
axB6.set_xlim([-0.05,0.2])
axB6.tick_params(axis='both', which='major', labelsize=8)
axB6.text(-0.35,0.93,'D',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axB6.transAxes)
hB3,lB3=axB6.get_legend_handles_labels()
hB3.remove(hB3[0])
lB3.remove(lB3[0])


axC6=sp.bivariable_std(collectr_Ispindles_duration,collectr_spindles_duration,axC6,stdx_zero=stdSpindlesb,meany_zero=np.mean(n_spindles_white),stdy_zero=np.std(n_spindles_white),marker='s',cmap=plt.cm.Reds ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_duration,meanzero_type='line_inf')
axC6.set_ylabel('$N_{SP}$/min.',fontsize=8,labelpad=0)
axC6.set_ylim([4,10])
axC6.set_xlim([-0.01,0.1])
axC6.tick_params(axis='both', which='major', labelsize=8)
axC6.text(-0.4,0.93,'E',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axC6.transAxes)

axD6=sp.bivariable_std(collectr_pso_duration,collectr_pcp_duration,axD6,meanx_zero=np.mean(P_SO_white),stdx_zero=np.std(P_SO_white),meany_zero=np.mean(P_PCP_white),stdy_zero=np.std(P_PCP_white),marker='s',cmap=plt.cm.Reds ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_duration,meanzero_type='line_inf')

axD6.set_ylabel('$P(C|SP)$',fontsize=8,labelpad=0)

axD6.set_ylim([0.1,0.5])
axD6.set_xlim([0.1,0.4])
axD6.plot([0,1],[0,1],color='gray')
axD6.set_xticks([0.1,0.2,0.3,0.4])
axD6.set_yticks([0.1,0.2,0.3,0.4,0.5])
axD6.tick_params(axis='both', which='major', labelsize=8)
axD6.text(-0.36,0.93,'F',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axD6.transAxes)


ax11r.set_axis_off()
ax11r.legend(hB1,lB1,fontsize=8,ncol=7,loc='upper left',bbox_to_anchor=(-0.01,0.1,1,1),columnspacing=3.0)

ax22r.set_axis_off()
legfig2=ax22r.legend(hB2,lB2,fontsize=8,ncol=7,loc='upper left',bbox_to_anchor=(-0.01,-.3,1,1),columnspacing=2.3)

ax33r.set_axis_off()
ax33r.legend(hB3,lB3,fontsize=8,ncol=7,loc='upper left',bbox_to_anchor=(-0.01,0.1,1,1),columnspacing=4.38)

fig2.savefig('./output/S2_Fig.eps',dpi=300,bbox_inches='tight',bbox_extra_artists=[legfig2])