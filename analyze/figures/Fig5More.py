#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 1 03:39:23 2021

@author: felipe
"""
import sys, os
sys.path.append(os.path.abspath(os.path.join('..', 'plotting')))
from cycler import cycler
import numpy as np
import scipy.stats as stats
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import statistic_plots as sp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import rc
rc('text', usetex=True)
timeSim=15.166
data_length=90999
colors=['#df009f','#008fdf','b','r','#ff7f0e','#666666','#f37782','#1f77b4']
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
RMSarray_white=fileSO['RMSarray_white']
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
stdSOb=np.std(np.sum(sWavelet[:,freqone:freqdot5],axis=1))
stdSpindlesb=np.std(np.sum(sWavelet[:,freqsixteen:freqeleven+1],axis=1))

#Load data
file_phase=np.load('../collectedData/CollectResults_phaseLongNew.npz')
collect_Iso_phase=file_phase['collect_Iso']
collect_Ispindles_phase=file_phase['collect_Ispindles']
collectr_Iso_phase=file_phase['collectr_Iso']
collectr_Ispindles_phase=file_phase['collectr_Ispindles']
collect_pc_phase=file_phase['collect_pc']
collect_pcp_phase=file_phase['collect_pcp']
collect_pso_phase=file_phase['collect_pso']
collect_psp_phase=file_phase['collect_psp']
collect_coincidences_phase=file_phase['collect_coincidences']
collect_so_phase=file_phase['collect_so']
collect_spindles_phase=file_phase['collect_spindles']
collectr_pc_phase=file_phase['collectr_pc']
collectr_pcp_phase=file_phase['collectr_pcp']
collectr_pso_phase=file_phase['collectr_pso']
collectr_psp_phase=file_phase['collectr_psp']
collectr_coincidences_phase=file_phase['collectr_coincidences']
collectr_so_phase=file_phase['collectr_so']
collectr_spindles_phase=file_phase['collectr_spindles']
y_axis_phase=[0,45,90]
labels_phase=['0','45','90']
labels_stimphase=['STIM-CL 0','STIM-CL 45','STIM-CL 90']

file_frequency=np.load('../collectedData/CollectResults_frequencyLong.npz')
collect_Iso_frequency=file_frequency['collect_Iso']
collect_Ispindles_frequency=file_frequency['collect_Ispindles']
collectr_Iso_frequency=file_frequency['collectr_Iso']
collectr_Ispindles_frequency=file_frequency['collectr_Ispindles']
collect_pc_frequency=file_frequency['collect_pc']
collect_pcp_frequency=file_frequency['collect_pcp']
collect_pso_frequency=file_frequency['collect_pso']
collect_psp_frequency=file_frequency['collect_psp']
collect_coincidences_frequency=file_frequency['collect_coincidences']
collect_so_frequency=file_frequency['collect_so']
collect_spindles_frequency=file_frequency['collect_spindles']
collectr_pc_frequency=file_frequency['collectr_pc']
collectr_pcp_frequency=file_frequency['collectr_pcp']
collectr_pso_frequency=file_frequency['collectr_pso']
collectr_psp_frequency=file_frequency['collectr_psp']
collectr_coincidences_frequency=file_frequency['collectr_coincidences']
collectr_so_frequency=file_frequency['collectr_so']
collectr_spindles_frequency=file_frequency['collectr_spindles']
y_axis_frequency=[0.5,0.65,0.85,1.05,1.25]
labels_frequency=['0.50','0.65','0.85','1.05','1.25']

#Figure Layout
fig=plt.figure(figsize=(5.2,1.7))
gs=gridspec.GridSpec(2,3,figure=fig,width_ratios=[0.5,0.5,0.5],height_ratios=[1,0.3])
gs.update(wspace=0.61, hspace=0.1) # set the spacing between axes.
axC=fig.add_subplot(gs[0,0])
axD=fig.add_subplot(gs[0,1])
axF=fig.add_subplot(gs[0,2])

#Plots
axC=sp.bivariable_std(collect_Iso_phase,collect_so_phase,axC,stdx_zero=stdSOb,meany_zero=np.mean(n_SO_white),stdy_zero=np.std(n_SO_white),marker='o',cmap=plt.cm.Dark2 ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_phase,meanzero_type='line_inf',colormap='direct')
axC=sp.bivariable_std(collectr_Iso_frequency[:,2:3],collectr_so_frequency[:,2:3],axC,stdx_zero=stdSOb,meany_zero=np.mean(n_SO_white),stdy_zero=np.std(n_SO_white),marker='|',cmap=colors[1],flagstdx=True,flagstdy=True,relative_diference=False,labels=['STIM-R'],meanzero_type='None')
axC=sp.bivariable_std(collect_Iso_frequency[:,2:3],collect_so_frequency[:,2:3],axC,stdx_zero=stdSOb,meany_zero=np.mean(n_SO_white),stdy_zero=np.std(n_SO_white),marker='_',cmap=colors[0] ,flagstdx=True,flagstdy=True,relative_diference=False,labels=['STIM-P'],meanzero_type='None')
axC.set_ylabel('$N_{SO}$/min.',fontsize=8)
axC.set_xlabel('$I^{(SO)}$ (a. u.)',fontsize=8)
axC.set_xlim([0.1,0.3])
axC.set_ylim([14,20])
axC.set_yticks([14,16,18,20])
axC.text(-0.4,1.01,'A',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axC.transAxes)


axD=sp.bivariable_std(collect_Ispindles_phase,collect_spindles_phase,axD,stdx_zero=stdSpindlesb,meany_zero=np.mean(n_spindles_white),stdy_zero=np.std(n_spindles_white),marker='o',cmap=plt.cm.Dark2 ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_phase,meanzero_type='line_inf',colormap='direct')
axD=sp.bivariable_std(collectr_Ispindles_frequency[:,2:3],collectr_spindles_frequency[:,2:3],axD,stdx_zero=stdSpindlesb,meany_zero=np.mean(n_spindles_white),stdy_zero=np.std(n_spindles_white),marker='|',cmap=colors[1] ,flagstdx=True,flagstdy=True,relative_diference=False,labels=['STIM-R'],meanzero_type='None')
axD=sp.bivariable_std(collect_Ispindles_frequency[:,2:3],collect_spindles_frequency[:,2:3],axD,stdx_zero=stdSpindlesb,meany_zero=np.mean(n_spindles_white),stdy_zero=np.std(n_spindles_white),marker='_',cmap=colors[0] ,flagstdx=True,flagstdy=True,relative_diference=False,labels=['STIM-P'],meanzero_type='None')
axD.set_ylabel('$N_{SP}$/min.',fontsize=8)
axD.set_xlabel('$I^{(SP)}$ (a. u.)')
axD.set_xlim([0.05,0.1])
axD.set_ylim([7,10])
axD.set_yticks([7,8,9,10])
axD.text(-0.4,1.01,'B',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axD.transAxes)


axF=sp.bivariable_std(collect_pso_phase,collect_pcp_phase,axF,meanx_zero=np.mean(P_SO_white),stdx_zero=np.std(P_SO_white),meany_zero=np.mean(P_PCP_white),stdy_zero=np.std(P_PCP_white),marker='o',cmap=plt.cm.Dark2,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_stimphase,meanzero_type='line_inf',colormap='direct')
axF=sp.bivariable_std(collectr_pso_frequency[:,2:3],collectr_pcp_frequency[:,2:3],axF,meanx_zero=np.mean(P_SO_white),stdx_zero=np.std(P_SO_white),meany_zero=np.mean(P_PCP_white),stdy_zero=np.std(P_PCP_white),marker='|',cmap=colors[1],flagstdx=True,flagstdy=True,relative_diference=False,labels=['STIM-R'],meanzero_type='None')
axF=sp.bivariable_std(collect_pso_frequency[:,2:3],collect_pcp_frequency[:,2:3],axF,meanx_zero=np.mean(P_SO_white),stdx_zero=np.std(P_SO_white),meany_zero=np.mean(P_PCP_white),stdy_zero=np.std(P_PCP_white),marker='_',cmap=colors[0] ,flagstdx=True,flagstdy=True,relative_diference=False,labels=['STIM-P'],meanzero_type='None')

axF.set_xlabel('$P(SO)$',fontsize=8)
axF.set_ylabel('$P(C|SP)$',fontsize=8)
axF.plot([0,1],[0,1],'gray')
axF.set_xlim([0.2,0.35])
axF.set_ylim([0.3,0.5])
axF.set_xticks([0.2,0.27,0.35])
axF.set_yticks([0.3,0.37,0.45])
axF.text(-0.48,1.01,'C',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axF.transAxes)
hF,lF=axF.get_legend_handles_labels()
axleg=fig.add_subplot(gs[1,:])
axleg.set_axis_off()
legF=axleg.legend(hF,lF,fontsize=8,ncol=3,loc='upper left',bbox_to_anchor=(0,-0.2,0.1,0.1),columnspacing=6)

#Save figure
fig.savefig('./output/Fig5New.eps',dpi=300,bbox_inches='tight',bbox_extra_artist=[legF])
