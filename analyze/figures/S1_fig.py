#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  9 04:05:11 2021

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
import CollectData.statistic_plots as sp
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

#Load stimulation data
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



#Figure layout
fig1=plt.figure(figsize=(5.2,3.5))
gs1=gridspec.GridSpec(4,3,figure=fig1,width_ratios=[1,1,1],height_ratios=[1,0.12,1,0.12])
gs1.update(wspace=0.5, hspace=0.4) # set the spacing between axes.
axB7=fig1.add_subplot(gs1[0,0])
axC7=fig1.add_subplot(gs1[0,1])
axD7=fig1.add_subplot(gs1[0,2])
axB8=fig1.add_subplot(gs1[2,0])
axC8=fig1.add_subplot(gs1[2,1])
axD8=fig1.add_subplot(gs1[2,2])
ax11=fig1.add_subplot(gs1[1,:])
ax22=fig1.add_subplot(gs1[3,:])

axB7=sp.bivariable_std(collect_Iso_frequency,collect_so_frequency,axB7,stdx_zero=stdSOb,meany_zero=np.mean(n_SO_white),stdy_zero=np.std(n_SO_white),marker='P',cmap=plt.cm.Blues ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_frequency,meanzero_type='line_inf')
axB7.set_ylabel('$N_{SO}$/min.',fontsize=8,labelpad=0)
axB7.set_ylim([12,22])
axB7.set_xlim([0.1,0.3])

# axB7.set_xticks([0.3,0.4,0.5,0.6,1])
# axB7.set_xticklabels([0.3,0.4,0.5,0.6,1],fontsize=8)
axB7.text(-0.36,0.9,'A',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axB7.transAxes)

hB,lB=axB7.get_legend_handles_labels()
hB7,lB7=axB7.get_legend_handles_labels()
hB7.remove(hB7[0])
lB7.remove(lB7[0])
for h in hB7:
    hB.append(h)
for l in lB7:
    lB.append(l)

axC7=sp.bivariable_std(collect_Ispindles_frequency,collect_spindles_frequency,axC7,stdx_zero=stdSpindlesb,meany_zero=np.mean(n_spindles_white),stdy_zero=np.std(n_spindles_white),marker='P',cmap=plt.cm.Blues ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_frequency,meanzero_type='line_inf')
axC7.set_ylim([8,10])
axC7.set_xlim([0.05,0.08])
axC7.set_ylabel('$N_{SP}$/min.',fontsize=8,labelpad=0)
# axC7.set_xticks([0.06,0.07,0.08,0.09,0.1,0.12])
# axC7.set_xticklabels([0.06,0.07,0.08,0.09,0.10,0.12],fontsize=8)
# axC7.legend(loc='upper left',fontsize=8)
axC7.text(-0.37,0.9,'B',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axC7.transAxes)

axD7=sp.bivariable_std(collect_pso_frequency,collect_pcp_frequency,axD7,meanx_zero=np.mean(P_SO_white),stdx_zero=np.std(P_SO_white),meany_zero=np.mean(P_PCP_white),stdy_zero=np.std(P_PCP_white),marker='P',cmap=plt.cm.Blues ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_frequency,meanzero_type='line_inf')
# axC1.set_title('Spindles')
axD7.set_ylabel('$P(C|SP)$',fontsize=8,labelpad=0)
# axC1.set_xscale('log')
axD7.set_ylim([0.25,0.55])
axD7.set_xlim([0.2,0.35])
axD7.plot([0,1],[0,1],color='gray')
axD7.set_yticks([0.25,0.40,0.55])
axD7.set_xticks([0.2,0.27,0.35])
axD7.text(-0.44,0.9,'C',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axD7.transAxes)



axB8=sp.bivariable_std(collectr_Iso_frequency,collectr_so_frequency,axB8,stdx_zero=stdSOb,meany_zero=np.mean(n_SO_white),stdy_zero=np.std(n_SO_white),marker='P',cmap=plt.cm.Purples ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_frequency,meanzero_type='line_inf')
axB8.set_ylabel('$N_{SO}$/min.',fontsize=8,labelpad=0)
axB8.set_xlabel('$I^{(SO)}$ (a. u.)',fontsize=8)
# axB8.set_xscale('log')
axB8.set_ylim([12,22])
axB8.set_xlim([0.1,0.3])
# axB8.set_xticks([0.3,0.4,0.5,0.6,1])
# axB8.set_xticklabels([0.3,0.4,0.5,0.6,1],fontsize=8)
axB8.text(-0.36,0.9,'D',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axB8.transAxes)

hB,lB=axB8.get_legend_handles_labels()
hB8,lB8=axB8.get_legend_handles_labels()
hB8.remove(hB8[0])
lB8.remove(lB8[0])
for h in hB8:
    hB.append(h)
for l in lB8:
    lB.append(l)

axC8=sp.bivariable_std(collectr_Ispindles_frequency,collectr_spindles_frequency,axC8,stdx_zero=stdSpindlesb,meany_zero=np.mean(n_spindles_white),stdy_zero=np.std(n_spindles_white),marker='P',cmap=plt.cm.Purples ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_frequency,meanzero_type='line_inf')
axC8.set_xlabel('$I^{(SP)}$ (a. u.)',fontsize=8,labelpad=0)
axC8.set_ylabel('$N_{SP}$/min.',fontsize=8,labelpad=0)
# axC8.set_xscale('log')
axC8.set_ylim([8,10])
axC8.set_xlim([0.05,0.08])
# axC8.set_xticks([0.06,0.07,0.08,0.09,0.1,0.12])
# axC8.set_xticklabels([0.06,0.07,0.08,0.09,0.10,0.12],fontsize=8)
# axC8.legend(loc='upper left',fontsize=8)
axC8.text(-0.37,0.9,'E',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axC8.transAxes)


axD8=sp.bivariable_std(collectr_pso_frequency,collectr_pcp_frequency,axD8,meanx_zero=np.mean(P_SO_white),stdx_zero=np.std(P_SO_white),meany_zero=np.mean(P_PCP_white),stdy_zero=np.std(P_PCP_white),marker='P',cmap=plt.cm.Purples ,flagstdx=True,flagstdy=True,relative_diference=False,labels=labels_frequency,meanzero_type='line_inf')
# axC1.set_title('Spindles')
axD8.set_ylabel('$P(C|SP)$',fontsize=8,labelpad=0)
# axC1.set_xscale('log')
axD8.set_ylim([0.25,0.55])
axD8.set_xlim([0.2,0.35])
axD8.plot([0,1],[0,1],color='gray')
axD8.set_yticks([0.25,0.40,0.55])
axD8.set_xticks([0.2,0.27,0.35])
axD8.set_xlabel('$P(SO)$',fontsize=8)
axD8.text(-0.44,0.9,'F',fontsize=12,fontweight=1000,verticalalignment='bottom',transform=axD8.transAxes)

ax11.set_axis_off()
legfig1=ax11.legend(hB7,lB7,fontsize=8,ncol=7,loc='upper left',bbox_to_anchor=(0,1.2,0.1,0.1),columnspacing=2.8)

ax22.set_axis_off()
leg2=ax22.legend(hB8,lB8,fontsize=8,ncol=7,loc='upper left',bbox_to_anchor=(0,-0.6,0.1,0.1),columnspacing=2.8)
fig1.savefig('./output/S1_fig.eps',dpi=300, bbox_extra_artists=[leg2])
