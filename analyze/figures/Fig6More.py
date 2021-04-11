#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 1 04:39:50 2021

@author: felipe
"""
import sys, os
sys.path.append(os.path.abspath(os.path.join('..', 'spectrum')))
sys.path.append(os.path.abspath(os.path.join('..', 'eventsDetection')))
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
fileSO=np.load('../collectedData/SWSbaselinesSOlongMore.npz')
SO_centers_white=fileSO['SO_centers_white']
SOarray_white=fileSO['SOarray_white']
RMSarray_white=fileSO['RMSarray_white']
c_white=int(fileSO['c_white'])
m_white=int(fileSO['m_white'])
tSO=np.arange(-2,2,0.01)

#Load data stimulation

file1=np.load('../collectedData/coincidenceSpindles/SWSlong-seed-5-decreaseRamp-F0.85-A3.46e+01-D1.00e-01-sdA-OFF-sdF-OFF-pF0.02-pR2.00-Phi.txtSpindlesCoincidenceLong.npz')
# file2=np.load('../../tables/coincidenceSpindles/SWSlong-seed-15-phase-0.00-decreaseRamp-F0.85-A3.46e+01-D1.00e-01-sdA-OFF-sdF-OFF-Phi.txtSpindlesCoincidenceLong.npz')
# file3=np.load('../../tables/coincidenceSpindles/SWSlong-seed-15-phase-0.79-decreaseRamp-F0.85-A3.46e+01-D1.00e-01-sdA-OFF-sdF-OFF-Phi.txtSpindlesCoincidenceLong.npz')
# file4=np.load('../../tables/coincidenceSpindles/SWSlong-seed-15-phase-1.57-decreaseRamp-F0.85-A3.46e+01-D1.00e-01-sdA-OFF-sdF-OFF-Phi.txtSpindlesCoincidenceLong.npz')
file2=np.load('/media/felipe/TOSHIBA2T/Phase-LongNew/tables/tables1/coincidenceSpindles/SWSLong-seed-5-phase-0.00-decreaseRamp-F0.85-A3.46e+01-D1.00e-01-sdA-OFF-sdF-ON-Phi.txtSpindlesCoincidenceLong.npz')
file3=np.load('/media/felipe/TOSHIBA2T/Phase-LongNew/tables/tables1/coincidenceSpindles/SWSLong-seed-5-phase-0.79-decreaseRamp-F0.85-A3.46e+01-D1.00e-01-sdA-OFF-sdF-ON-Phi.txtSpindlesCoincidenceLong.npz')
file4=np.load('/media/felipe/TOSHIBA2T/Phase-LongNew/tables/tables1/coincidenceSpindles/SWSLong-seed-5-phase-1.57-decreaseRamp-F0.85-A3.46e+01-D1.00e-01-sdA-OFF-sdF-ON-Phi.txtSpindlesCoincidenceLong.npz')


n_spindles=file1['n_spindles']/timeSim
n_Delta=file1['n_Delta']/timeSim
n_SO=file1['n_SO']/timeSim
n_SOp=file1['n_SOp']/timeSim
n_spindlesp=file1['n_spindlesp']/timeSim
n_Deltap=file1['n_Deltap']/timeSim
coincident_centers=file1['coincident_centers']
coincident_centersp=file1['coincident_centersp']
SO_Dpeak_phase=file1['SO_Dpeak_phase']
SO_Dpeak_phasep=file1['SO_Dpeak_phasep']
SO_Dpeak_time=file1['SO_Dpeak_time']
SO_Dpeak_timep=file1['SO_Dpeak_timep']
SpindlesSOcount=file1['SpindlesSOcount']
SpindlesSOcountp=file1['SpindlesSOcountp']
Spindles_centers_phase=file1['Spindles_centers_phase']
Spindles_centers_phasep=file1['Spindles_centers_phasep']
Spindles_coincidences_phase=file1['Spindles_coincidences_phase']
Spindles_coincidences_time=file1['Spindles_coincidences_time']
Spindles_coincidences_phasep=file1['Spindles_coincidences_phasep']
Spindles_coincidences_timep=file1['Spindles_coincidences_timep']
Spindles_centers=file1['Spindles_centers']
Spindles_centersp=file1['Spindles_centersp']
SOarray=file1['SOarray']
RMSarray=file1['RMSarray']
SOarrayp=file1['SOarrayp']
RMSarrayp=file1['RMSarrayp']
SOextracted=file1['SOextracted']
SOextractedp=file1['SOextractedp']

n_spindlesPh=file2['n_spindles']/timeSim
n_DeltaPh=file2['n_Delta']/timeSim
n_SOPh=file2['n_SO']/timeSim
n_SOpPh=file2['n_SOp']/timeSim
n_spindlespPh=file2['n_spindlesp']/timeSim
n_DeltapPh=file2['n_Deltap']/timeSim
coincident_centersPh=file2['coincident_centers']
coincident_centerspPh=file2['coincident_centersp']
SO_Dpeak_phasePh=file2['SO_Dpeak_phase']
SO_Dpeak_phasepPh=file2['SO_Dpeak_phasep']
SO_Dpeak_timePh=file2['SO_Dpeak_time']
SO_Dpeak_timepPh=file2['SO_Dpeak_timep']
SpindlesSOcountPh=file2['SpindlesSOcount']
SpindlesSOcountpPh=file2['SpindlesSOcountp']
Spindles_centers_phasePh=file2['Spindles_centers_phase']
Spindles_centers_phasepPh=file2['Spindles_centers_phasep']
Spindles_coincidences_phasePh=file2['Spindles_coincidences_phase']
Spindles_coincidences_timePh=file2['Spindles_coincidences_time']
Spindles_coincidences_phasepPh=file2['Spindles_coincidences_phasep']
Spindles_coincidences_timepPh=file2['Spindles_coincidences_timep']
Spindles_centersPh=file2['Spindles_centers']
Spindles_centerspPh=file2['Spindles_centersp']
SOarrayPh=file2['SOarray']
SOarraypPh=file2['SOarrayp']
RMSarrayPh=file2['RMSarray']
SOextractedPh=file2['SOextracted']
SOextractedpPh=file2['SOextractedp']


n_spindlesPh1=file3['n_spindles']/timeSim
n_DeltaPh1=file3['n_Delta']/timeSim
n_SOPh1=file3['n_SO']/timeSim
n_SOpPh1=file3['n_SOp']/timeSim
n_spindlespPh1=file3['n_spindlesp']/timeSim
n_DeltapPh1=file3['n_Deltap']/timeSim
coincident_centersPh1=file3['coincident_centers']
coincident_centerspPh1=file3['coincident_centersp']
SO_Dpeak_phasePh1=file3['SO_Dpeak_phase']
SO_Dpeak_phasepPh1=file3['SO_Dpeak_phasep']
SO_Dpeak_timePh1=file3['SO_Dpeak_time']
SO_Dpeak_timepPh1=file3['SO_Dpeak_timep']
SpindlesSOcountPh1=file3['SpindlesSOcount']
SpindlesSOcountpPh1=file3['SpindlesSOcountp']
Spindles_centers_phasePh1=file3['Spindles_centers_phase']
Spindles_centers_phasepPh1=file3['Spindles_centers_phasep']
Spindles_coincidences_phasePh1=file3['Spindles_coincidences_phase']
Spindles_coincidences_timePh1=file3['Spindles_coincidences_time']
Spindles_coincidences_phasepPh1=file3['Spindles_coincidences_phasep']
Spindles_coincidences_timepPh1=file3['Spindles_coincidences_timep']
Spindles_centersPh1=file3['Spindles_centers']
Spindles_centerspPh1=file3['Spindles_centersp']
SOarrayPh1=file3['SOarray']
SOarraypPh1=file3['SOarrayp']
RMSarrayPh1=file3['RMSarray']
SOextractedPh1=file3['SOextracted']
SOextractedpPh1=file3['SOextractedp']


n_spindlesPh2=file4['n_spindles']/timeSim
n_DeltaPh2=file4['n_Delta']/timeSim
n_SOPh2=file4['n_SO']/timeSim
n_SOpPh2=file4['n_SOp']/timeSim
n_spindlespPh2=file4['n_spindlesp']/timeSim
n_DeltapPh2=file4['n_Deltap']/timeSim
coincident_centersPh2=file4['coincident_centers']
coincident_centerspPh2=file4['coincident_centersp']
SO_Dpeak_phasePh2=file4['SO_Dpeak_phase']
SO_Dpeak_phasepPh2=file4['SO_Dpeak_phasep']
SO_Dpeak_timePh2=file4['SO_Dpeak_time']
SO_Dpeak_timepPh2=file4['SO_Dpeak_timep']
SpindlesSOcountPh2=file4['SpindlesSOcount']
SpindlesSOcountpPh2=file4['SpindlesSOcountp']
Spindles_centers_phasePh2=file4['Spindles_centers_phase']
Spindles_centers_phasepPh2=file4['Spindles_centers_phasep']
Spindles_coincidences_phasePh2=file4['Spindles_coincidences_phase']
Spindles_coincidences_timePh2=file4['Spindles_coincidences_time']
Spindles_coincidences_phasepPh2=file4['Spindles_coincidences_phasep']
Spindles_coincidences_timepPh2=file4['Spindles_coincidences_timep']
Spindles_centersPh2=file4['Spindles_centers']
Spindles_centerspPh2=file4['Spindles_centersp']
SOarrayPh2=file4['SOarray']
SOarraypPh2=file4['SOarrayp']
RMSarrayPh2=file4['RMSarray']
SOextractedPh2=file4['SOextracted']
SOextractedpPh2=file4['SOextractedp']

#Figure layout
fig6=plt.figure(figsize=(5.2,5.2))
gs6=gridspec.GridSpec(1,1,figure=fig6,width_ratios=[1],height_ratios=[1])
axE=fig6.add_subplot(gs6[0])

Spindles_coincidences_time_baseline=Spindles_coincenters_time_white[0:SpindlesSOcount_whiteall]
Spindles_coincidences_time_all=Spindles_coincidences_time[0,0:SpindlesSOcount[0]]
Spindles_coincidences_timer_all=Spindles_coincidences_timep[0,0:SpindlesSOcountp[0]]
Spindles_coincidences_timep_all=Spindles_coincidences_timePh[0,0:SpindlesSOcountPh[0]]
Spindles_coincidences_timep1_all=Spindles_coincidences_timePh1[0,0:SpindlesSOcountPh1[0]]
Spindles_coincidences_timep2_all=Spindles_coincidences_timePh2[0,0:SpindlesSOcountPh2[0]]

SOarrayC=SOarray[0,0:SOextracted[0],:]
SOarrayR=SOarrayp[0,0:SOextractedp[0],:]
SOarrayP=SOarrayPh[0,0:SOextractedPh[0],:]
SOarrayP1=SOarrayPh1[0,0:SOextractedPh1[0],:]
SOarrayP2=SOarrayPh2[0,0:SOextractedPh2[0],:]
RMSarrayC=RMSarray[0,0:SOextracted[0],:]
RMSarrayR=RMSarrayp[0,0:SOextractedp[0],:]
RMSarrayP=RMSarrayPh[0,0:SOextractedPh[0],:]
RMSarrayP1=RMSarrayPh1[0,0:SOextractedPh1[0],:]
RMSarrayP2=RMSarrayPh2[0,0:SOextractedPh2[0],:]
for seed in range(1,5):
    Spindles_coincidences_time_all=np.hstack((Spindles_coincidences_time_all,Spindles_coincidences_time[seed,0:SpindlesSOcount[seed]]))
    Spindles_coincidences_timer_all=np.hstack((Spindles_coincidences_timer_all,Spindles_coincidences_timep[seed,0:SpindlesSOcountp[seed]]))
    Spindles_coincidences_timep_all=np.hstack((Spindles_coincidences_timep_all,Spindles_coincidences_timePh[seed,0:SpindlesSOcountPh[seed]]))
    Spindles_coincidences_timep1_all=np.hstack((Spindles_coincidences_timep1_all,Spindles_coincidences_timePh1[seed,0:SpindlesSOcountPh1[seed]]))
    Spindles_coincidences_timep2_all=np.hstack((Spindles_coincidences_timep2_all,Spindles_coincidences_timePh2[seed,0:SpindlesSOcountPh2[seed]]))

    SOarrayC=np.vstack((SOarrayC,SOarray[seed,0:SOextracted[seed],:]))
    SOarrayR=np.vstack((SOarrayR,SOarrayp[seed,0:SOextractedp[seed],:]))
    SOarrayP=np.vstack((SOarrayP,SOarrayPh[seed,0:SOextractedPh[seed],:]))
    SOarrayP1=np.vstack((SOarrayP1,SOarrayPh1[seed,0:SOextractedPh1[seed],:]))
    SOarrayP2=np.vstack((SOarrayP2,SOarrayPh1[seed,0:SOextractedPh2[seed],:]))
    RMSarrayC=np.vstack((RMSarrayC,RMSarray[seed,0:SOextracted[seed],:]))
    RMSarrayR=np.vstack((RMSarrayR,RMSarrayp[seed,0:SOextractedp[seed],:]))
    RMSarrayP=np.vstack((RMSarrayP,RMSarrayPh[seed,0:SOextractedPh[seed],:]))
    RMSarrayP1=np.vstack((RMSarrayP1,RMSarrayPh1[seed,0:SOextractedPh1[seed],:]))
    RMSarrayP2=np.vstack((RMSarrayP2,RMSarrayPh2[seed,0:SOextractedPh2[seed],:]))
#Mean of SO
meanSOarrayC=np.mean(SOarrayC[0:np.sum(SOextracted),:],axis=0)
meanSOarrayR=np.mean(SOarrayR[0:np.sum(SOextractedp),:],axis=0)
meanSOarrayP=np.mean(SOarrayP[0:np.sum(SOextractedPh),:],axis=0)
meanSOarrayP1=np.mean(SOarrayP1[0:np.sum(SOextractedPh1),:],axis=0)
meanSOarrayP2=np.mean(SOarrayP2[0:np.sum(SOextractedPh2),:],axis=0)
meanSOarrayB=np.mean(SOarray_white[0:c_white,:],axis=0)


#Sort by amplitude
baseline_order=np.flip(np.argsort(np.max(np.abs(SOarray_white),axis=1)))[0:100]
c_order=np.flip(np.argsort(np.max(np.abs(SOarrayC),axis=1)))[0:100]
r_order=np.flip(np.argsort(np.max(np.abs(SOarrayR),axis=1)))[0:100]
p_order=np.flip(np.argsort(np.max(np.abs(SOarrayP),axis=1)))[0:100]
p1_order=np.flip(np.argsort(np.max(np.abs(SOarrayP1),axis=1)))[0:100]
p2_order=np.flip(np.argsort(np.max(np.abs(SOarrayP2),axis=1)))[0:100]

tBC_SO,pBC_SO=stats.ttest_ind(SOarray_white[baseline_order,:],SOarrayC[c_order,:])
tBR_SO,pBR_SO=stats.ttest_ind(SOarray_white[baseline_order,:],SOarrayR[r_order,:])
tBP_SO,pBP_SO=stats.ttest_ind(SOarray_white[baseline_order,:],SOarrayP[p_order,:])
tBP1_SO,pBP1_SO=stats.ttest_ind(SOarray_white[baseline_order,:],SOarrayP1[p1_order,:])
tBP2_SO,pBP2_SO=stats.ttest_ind(SOarray_white[baseline_order,:],SOarrayP2[p2_order,:])
tCP_SO,pCP_SO=stats.ttest_ind(SOarrayC[c_order,:],SOarrayP[p_order,:])
tCR_SO,pCR_SO=stats.ttest_ind(SOarrayC[c_order,:],SOarrayR[r_order,:])
tRP_SO,pRP_SO=stats.ttest_ind(SOarrayR[r_order,:],SOarrayP[p_order,:])

tBC_SP,pBC_SP=stats.ttest_ind(RMSarray_white[baseline_order,:],RMSarrayC[c_order,:])
tBR_SP,pBR_SP=stats.ttest_ind(RMSarray_white[baseline_order,:],RMSarrayR[r_order,:])
tBR_SP,pBP_SP=stats.ttest_ind(RMSarray_white[baseline_order,:],RMSarrayP[p_order,:])
tBR_SP,pBP1_SP=stats.ttest_ind(RMSarray_white[baseline_order,:],RMSarrayP1[p_order,:])
tBR_SP,pBP2_SP=stats.ttest_ind(RMSarray_white[baseline_order,:],RMSarrayP2[p_order,:])
tCP_SP,pCP_SP=stats.ttest_ind(RMSarrayC[c_order,:],RMSarrayP[p_order,:])
tCR_SP,pCR_SP=stats.ttest_ind(RMSarrayC[c_order,:],RMSarrayR[r_order,:])
tRP_SP,pRP_SP=stats.ttest_ind(RMSarrayR[r_order,:],RMSarrayP[p_order,:])

pOnes=np.ones((400,))
pBC_SOMasked=np.ma.masked_where(pBC_SO>=0.01, pOnes)*56
pBR_SOMasked=np.ma.masked_where(pBR_SO>=0.01, pOnes)*57
pBP_SOMasked=np.ma.masked_where(pBP_SO>=0.01, pOnes)*58
pBP1_SOMasked=np.ma.masked_where(pBP1_SO>=0.01, pOnes)*59
pBP2_SOMasked=np.ma.masked_where(pBP2_SO>=0.01, pOnes)*60
# pCP_SOMasked=np.ma.masked_where(pCP_SO>=0.01, pOnes)*48
# pCR_SOMasked=np.ma.masked_where(pCR_SO>=0.01, pOnes)*49
# pRP_SOMasked=np.ma.masked_where(pRP_SO>=0.01, pOnes)*50


pOnes=np.ones((400,))
pBC_RMSMasked=np.ma.masked_where(pBC_SP>=0.01, pOnes)*0.61
pBR_RMSMasked=np.ma.masked_where(pBR_SP>=0.01, pOnes)*0.62
pBP_RMSMasked=np.ma.masked_where(pBP_SP>=0.01, pOnes)*0.63
pCP_RMSMasked=np.ma.masked_where(pCP_SP>=0.01, pOnes)*0.64
pCR_RMSMasked=np.ma.masked_where(pCR_SP>=0.01, pOnes)*0.65
pRP_RMSMasked=np.ma.masked_where(pRP_SP>=0.01, pOnes)*0.66

#Wilcoxon test
tRB,pRB=stats.wilcoxon(Spindles_coincidences_timer_all[0:77],Spindles_coincidences_time_baseline[0:77],zero_method='pratt')
tCB,pCB=stats.wilcoxon(Spindles_coincidences_time_all[0:77],Spindles_coincidences_time_baseline[0:77],zero_method='pratt')
tPB,pPB=stats.wilcoxon(Spindles_coincidences_timep_all[0:77],Spindles_coincidences_time_baseline[0:77],zero_method='pratt')
tPB1,pPB1=stats.wilcoxon(Spindles_coincidences_timep1_all[0:77],Spindles_coincidences_time_baseline[0:77],zero_method='pratt')
tPB2,pPB2=stats.wilcoxon(Spindles_coincidences_timep2_all[0:77],Spindles_coincidences_time_baseline[0:77],zero_method='pratt')

tCR,pCR=stats.wilcoxon(Spindles_coincidences_time_all[0:77],Spindles_coincidences_timer_all[0:77],zero_method='pratt')
tRP,pRP=stats.wilcoxon(Spindles_coincidences_timer_all[0:77],Spindles_coincidences_timep_all[0:77],zero_method='pratt')
tCP,pCP=stats.wilcoxon(Spindles_coincidences_time_all[0:77],Spindles_coincidences_timep_all[0:77],zero_method='pratt')
print('Spindles distribution of delays')
print(pCB,pRB,pPB,pPB1,pPB2,pCR,pRP,pCP)

tRSk,pRSk=stats.ks_2samp(Spindles_coincidences_timer_all,Spindles_coincidences_time_baseline)
tPSk,pPSk=stats.ks_2samp(Spindles_coincidences_time_all,Spindles_coincidences_time_baseline)
tRC1Sk,pC1Sk=stats.ks_2samp(Spindles_coincidences_timep_all,Spindles_coincidences_time_baseline)
tC2Sk,pC2Sk=stats.ks_2samp(Spindles_coincidences_timep1_all,Spindles_coincidences_time_baseline)
tPC3Sk,pC3Sk=stats.ks_2samp(Spindles_coincidences_timep2_all,Spindles_coincidences_time_baseline)

tPRk,pPRk=stats.ks_2samp(Spindles_coincidences_time_all,Spindles_coincidences_timer_all)
tC1Rk,pC1Rk=stats.ks_2samp(Spindles_coincidences_timep_all,Spindles_coincidences_timer_all)
tC2Rk,pC2Rk=stats.ks_2samp(Spindles_coincidences_timep1_all,Spindles_coincidences_timer_all)
tC3Rk,pC3Rk=stats.ks_2samp(Spindles_coincidences_timep2_all,Spindles_coincidences_timer_all)

tC1Pk,pC1Pk=stats.ks_2samp(Spindles_coincidences_timep_all,Spindles_coincidences_time_all)
tC2Pk,pC2Pk=stats.ks_2samp(Spindles_coincidences_timep1_all,Spindles_coincidences_time_all)
tC3Pk,pC3Pk=stats.ks_2samp(Spindles_coincidences_timep2_all,Spindles_coincidences_time_all)


tC2C1k,pC2C1k=stats.ks_2samp(Spindles_coincidences_timep1_all,Spindles_coincidences_timep_all)
tC3C1k,pC3C1k=stats.ks_2samp(Spindles_coincidences_timep2_all,Spindles_coincidences_timep_all)

tC3C2k,pC3C2k=stats.ks_2samp(Spindles_coincidences_timep2_all,Spindles_coincidences_timep1_all)

print('Spindles distribution kolmogorov')
print(pRSk,pPSk,pC1Sk,pC2Sk,pC3Sk)
print(pPRk, pC1Rk, pC2Rk, pC3Rk)
print(pC1Pk, pC2Pk, pC3Pk)
print(pC2C1k, pC3C1k)
print(pC3C2k)


closed_loop_times=[-1*Spindles_coincidences_timep_all,-1*Spindles_coincidences_timep1_all,-1*Spindles_coincidences_timep2_all]

colorsHist=[cm.Dark2(0),cm.Dark2(1),cm.Dark2(2)]
axE.hist(closed_loop_times,bins=np.arange(-1.55,1.6,0.1),histtype='bar',color=colorsHist ,density=True)
axE.hist(-1*Spindles_coincidences_time_baseline,bins=np.arange(-1.55,1.6,0.1),histtype='step',color='black',linewidth=1.5,density=True)
axE1 = axE.twinx()
axE1.plot(tSO,meanSOarrayP*1,color=cm.Dark2(0),label='STIM-CL 0')
axE1.plot(tSO,meanSOarrayP1*1,color=cm.Dark2(1),label='STIM-CL 45')
axE1.plot(tSO,meanSOarrayP2*1,color=cm.Dark2(2),label='STIM-CL 90')
axE1.plot(tSO,meanSOarrayB*1,'--',color='black',label='SHAM')
axE.set_ylabel('pdf',fontsize=8)
axE1.set_ylabel('SO amplitude ($s^{-1}$)',fontsize=8)
axE.set_xlabel('Delay from SO DOWN-peak (s)',fontsize=8)
axE.set_ylim([0,2.5])
axE.set_xlim([-1.6,1.6])
axE1.set_ylim([-2.5,1])
leg1=axE1.legend(ncol=4,loc='upper left',fontsize=8,bbox_to_anchor=(0,0,0,1),columnspacing=1.4)

fig6.savefig('./output/Fig6New.eps',dpi=300,bbox_inches='tight',bbox_extra_artist=[leg1])