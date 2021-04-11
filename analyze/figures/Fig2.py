
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot of Fig2.eps
Created on Sat Jul 18 14:51:56 2020
Modified on Mon Feb 1 13:14:12 2021
@author: felipe
@project: Selection of stimulus parameters for enhancing slow wave sleep events with a Neural-field theory thalamocortical computational model
"""
import sys, os
sys.path.append(os.path.abspath(os.path.join('..', 'spectrum')))
sys.path.append(os.path.abspath(os.path.join('..', 'eventsDetection')))
sys.path.append(os.path.abspath(os.path.join('..', 'plotting')))
import numpy as np
import scipy.signal as signal
import scipy.stats as stats
from scipy import ndimage
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
import pandas as pd
import Data
import spindles
import SO
import ERP
import Wavelets
import re
import h5py
import steadyFunctions

from matplotlib import rc
rc('text', usetex=True)

colorsa=['#df009f','#008fdf','b','r','#ff7f0e','#666666','#f37782','#1f77b4']
h1=1e-2
fs=100
#Load timeseries files
filename = '/media/felipe/TOSHIBA2T/Frequency-White/SWS-seed-1-halfTrapece-F0.85-A2.45e+01-D1.00e-01-sdA-OFF-sdF-OFF-pF0.02-pR2.00-Phi.txt'
filenameCL = '/media/felipe/TOSHIBA2T/Phase-White/SWS-seed-1-phase-0.79-halfTrapece-F0.85-A2.45e+01-D1.00e-01-sdA-OFF-sdF-OFF-pF0.02-pR2.00-Phi.txt'
filenamePoisson=filename.replace('sdF-OFF','sdF-Poisson')
filenamePhase='/media/felipe/TOSHIBA2T/Phase-White/SWS-seed-1-phase-0.79-halfTrapece-F0.85-A2.45e+01-D1.00e-01-sdA-OFF-sdF-OFF-pF0.02-pR2.00-Phase.txt'
indexPath=0		
for m in range(len(filename)-1,0,-1):
    if filename[m]=='/':
        indexPath=m
        break
np.random.seed(0)

#Baseline of the same seed
nameTokens=filename.split('-')
print(filename)
filenameBaseline='/media/felipe/TOSHIBA2T/Baselines-White/SWS-seed-'+nameTokens[3]+'-baseline-Phi.txt' 

#Load Phi, Stim and Markers
Phie1,Phiez,Phinz,Baselineez,Baselinen,Stim,marker,time=Data.loadAllData(filename, filenameBaseline)
Phieph1,Phiephz,Phinphz,Baselineephz,Baselinephn,Stimph,markerph,timeph=Data.loadAllData(filenameCL, filenameBaseline)
#Load baseline data
Ee,Ez,stdB=Data.loadStoredData1Ch(filenameBaseline,L=np.arange(0,256))
#Load aditional data
# Re,Rz,stdR=Data.loadStoredData1Ch(filenameCL,L=[512])
#Load phase
phase_info,phasez,stdphase=Data.loadStoredData1Ch(filenamePhase,L=[0,1,2],sep=',')
Phiep1,Phiepz,Phipnz,Baselineep,Baselinen,Stim,markers,time=Data.loadAllData(filenamePoisson, filenameBaseline)
phase_info=phase_info.to_numpy()
phase=phase_info[:,2]
other_online=phase_info[:,1]
delta_online=phase_info[:,0]
#Calculate z-score of baseline
Ee=Ee.to_numpy()
Ee=np.mean(Ee,axis=1)-np.mean(Ee)
print(filename[indexPath::])
print('Loaded data')
dataStim=np.nanmean(Stim,axis=1)

#%%
#SP Filter 
bSp,aSp=signal.cheby1(4,1e-6,[9/(fs/2), 16/(fs/2)],'bandpass')
filtered_S=signal.filtfilt(bSp,aSp,Baselineez)
stdSp=np.std(filtered_S)
meanRMSSp=spindles.meanRMS(Baselineez,fs=fs)
#SO Filter
bD,aD=signal.cheby1(4,1e-6,[0.5/(fs/2), 1.25/(fs/2)],'bandpass')
filtered_D=signal.filtfilt(bD,aD,Baselineez)
stdD=np.std(filtered_D)
#Theta filter
bT,aT=signal.cheby1(4,1e-6,[4.5/(fs/2), 8/(fs/2)],'bandpass')
filtered_T=signal.filtfilt(bT,aT,Baselineez)
stdT=np.std(filtered_T)

b85,a85=signal.cheby1(4,0.001,[0.83/(fs/2), 0.87/(fs/2)],'bandpass')
filtered_85=signal.filtfilt(b85,a85,Phieph1)

print('Baseline')
detectionSO,Sosegmentation, SOs,countSO,meanSOb=SO.SODetection(Ee,neg_threshold=-1e-6,p2p_threshold=0,SOp2p_threshold=0)
filtered, rms,detectionSpindles, countSpindles=spindles.detectSpindles(Baselineez,fs=fs,threshold=1.25*meanRMSSp)  
filtered, rms,detectionTheta, countTheta=spindles.detectBand(Baselineez,fs=fs,lowFrequency=4.5,highFrequency=8.0,threshold=1.2*stdT,tmin=0.5,gap_max=0.05,windowLength=0.5)
filtered, rms,detectionDelta, countDelta=spindles.detectBand(Baselineez,fs=fs,lowFrequency=0.5,highFrequency=1.25,threshold=1.2*stdD,tmin=1.0,gap_max=0.05,windowLength=0.5)

print('rhythmic')
detectionSOC,Sosegmentationc, SOs,countSO,mSo=SO.SODetection(Phie1,neg_threshold=-1e-6,p2p_threshold=0,SOp2p_threshold=1.25*meanSOb)
filtered, rms,detectionSpindlesC, countSpindlesc=spindles.detectSpindles(Phiez,fs=fs,threshold=1.25*meanRMSSp)  
filtered, rms,detectionThetaC, countTheta=spindles.detectBand(Phiez,fs=fs,lowFrequency=4.5,highFrequency=8.0,threshold=1.2*stdT,tmin=0.5,gap_max=0.05,windowLength=0.5)
filtered, rms,detectionDeltaC, countDelta=spindles.detectBand(Phiez,fs=fs,lowFrequency=0.5,highFrequency=1.25,threshold=1.2*stdD,tmin=1.0,gap_max=0.05,windowLength=0.5)
print('random')
detectionSOR,Sosegmentationr,SOs,countSO,mSo=SO.SODetection(Phiep1,neg_threshold=-1e-6,p2p_threshold=0,SOp2p_threshold=1.25*meanSOb)
filtered, rms,detectionSpindlesR, countSpindlesr=spindles.detectSpindles(Phiepz,fs=fs,threshold=1.25*meanRMSSp)  
filtered, rms,detectionThetaR, countTheta=spindles.detectBand(Phiepz,fs=fs,lowFrequency=4.5,highFrequency=8.0,threshold=1.2*stdT,tmin=0.5,gap_max=0.05,windowLength=0.5)
filtered, rms,detectionDeltaR, countDelta=spindles.detectBand(Phiepz,fs=fs,lowFrequency=0.5,highFrequency=1.25,threshold=1.2*stdD,tmin=1.0,gap_max=0.05,windowLength=0.5)

print('Phase')
detectionSOP,Sosegmentationp, SOsP,countSOP,mSo=SO.SODetection(Phieph1,neg_threshold=-1e-6,p2p_threshold=0,SOp2p_threshold=1.25*meanSOb)
filtered_sp, rms,detectionSpindlesP, countSpindlesp=spindles.detectSpindles(Phiephz,fs=fs,threshold=1.25*meanRMSSp)  
filtered, rms,detectionThetaP, countTheta=spindles.detectBand(Phiephz,fs=fs,lowFrequency=4.5,highFrequency=8.0,threshold=1.2*stdT,tmin=0.5,gap_max=0.05,windowLength=0.5)
filtered_d, rms,detectionDeltaP, countDelta=spindles.detectBand(Phiephz,fs=fs,lowFrequency=0.5,highFrequency=1.25,threshold=1.2*stdD,tmin=1.0,gap_max=0.05,windowLength=0.5)


print(countSpindles)

# print('Welch periodograms')
# #Frequency
# f,P=Data.welchSingle(Phiez)
# f,Prandom=Data.welchSingle(Phiepz)
# f,Pclosed=Data.welchSingle(Phiephz)
# f,Pstim=Data.welchSingle(dataStim)
# f,Pn=Data.welchSingle(Phinz)
# f,Pbaseline=Data.welchSingle(Baselineez)

fileEvents=np.load('../collectedData/EventsCharacteristicsBaseline.npz')
fileEventsPR=np.load('../collectedData/EventsCharacteristicsPR.npz')
durationsSOB=fileEvents['durationsSO']
durationsSOP=fileEventsPR['durationsSOP']
durationsSOR=fileEventsPR['durationsSOR']
durationsSPB=fileEvents['durationsSP']
durationsSPP=fileEventsPR['durationsSPP']
durationsSPR=fileEventsPR['durationsSOR']
p2pSOB=fileEvents['p2pSO']
p2pSOP=fileEventsPR['p2pSOP']
p2pSOR=fileEventsPR['p2pSOR']
p2pSPB=fileEvents['p2pSP']
p2pSPP=fileEventsPR['p2pSPP']
p2pSPR=fileEventsPR['p2pSPR']


print('Scalograms')
freqsWavelet,scales,Scalogram=Data.waveletSingle(Phiez,correctF=False)
#freqsWavelet,scales,ScalogramC=Data.waveletSingle(Phiez,correctF=True)
freqsWavelet,scales,ScalogramRandom=Data.waveletSingle(Phiepz,correctF=False)
#freqsWavelet,scales,ScalogramRandomC=Data.waveletSingle(Phiepz,correctF=True)
freqB,scalesB,ScalogramB=Data.waveletSingle(Baselineez,correctF=False)
# freqB,scalesB,ScalogramBC=Data.waveletSingle(Baselineez,correctF=True)
freqsPhase,scalesPhase,ScalogramPhase=Data.waveletSingle(Phiephz,correctF=False)
# freqsPhase,scalesPhase,ScalogramPhaseC=Data.waveletSingle(Phiephz,correctF=True)
#%%
prop_cycle = plt.rcParams['axes.prop_cycle']
colorsProp = prop_cycle.by_key()['color']
#Figure layout
fig1=plt.figure(figsize=(8.82,8.65))
gs1 = gridspec.GridSpec(12, 3, figure=fig1, 
                        height_ratios=[0.2,0.2,0.2,0.2,0.1,0.2,0.04,0.25,0.04,0.2,0.04,0.2],
                        width_ratios=[1.2,1,0.05],wspace=0.4,hspace=0.3)
#Time series
axA1= fig1.add_subplot(gs1[0,0])
axA2= fig1.add_subplot(gs1[1,0])
axA3= fig1.add_subplot(gs1[2,0])
axA4= fig1.add_subplot(gs1[3,0])
#Spectrograms
axB1= fig1.add_subplot(gs1[0,1])
axB2= fig1.add_subplot(gs1[1,1])
axB3= fig1.add_subplot(gs1[2,1])
axB4= fig1.add_subplot(gs1[3,1])
#Zoom
axC1= fig1.add_subplot(gs1[5,0])
axC2= fig1.add_subplot(gs1[7,0])
axC3= fig1.add_subplot(gs1[9,0])
#Shapes
#Spectrum 
# axD= fig1.add_subplot(gs1[5,1:3])
axE= fig1.add_subplot(gs1[5:10,1:3])

axB5=fig1.add_subplot(gs1[0:4,2])

# axD1= fig1.add_subplot(gs1[4,1:3])
# axD1.set_axis_off()

axC4=fig1.add_subplot(gs1[5:10,2])
axC4.set_axis_off()
axC5=fig1.add_subplot(gs1[6,0])
axC5.set_axis_off()
axC6=fig1.add_subplot(gs1[8,0])
axC6.set_axis_off()

#Events
axF1=fig1.add_subplot(gs1[10,:])
axF1.set_axis_off()

axF=fig1.add_subplot(gs1[11,0])
axG=fig1.add_subplot(gs1[11,1:3])


#Shape
A=np.array([10,17.32,17.32,17.32,15.897,15.885,14.883,14.142,16.330,14.142,20.186,30.001,12.248,12.248])
A=A*2
duration=0.1
fs=1000


indexesWavelet=[33,69,106,142,178,215]
freqs=freqsWavelet[indexesWavelet]
labels_appFreq=[0.5,1,2,4,8,16]
lenF,lenT=np.shape(Scalogram)
initT=1500
lenT=3000
fs=100
nticks=30

offset=0.125
spindlesMasked=np.ma.masked_where(detectionSpindles==0, detectionSpindles==1)*(-0.032)
SOMasked=np.ma.masked_where(detectionSO==0, detectionSO==1)*(0.032)
spindlesMaskedP=np.ma.masked_where(detectionSpindlesP==0, detectionSpindlesP==1)*(-0.032)
SOMaskedP=np.ma.masked_where(detectionSOP==0, detectionSOP==1)*(0.032)
spindlesMaskedC=np.ma.masked_where(detectionSpindlesC==0, detectionSpindlesC==1)*(-0.032)
SOMaskedC=np.ma.masked_where(detectionSOC==0, detectionSOC==1)*(0.032)
spindlesMaskedR=np.ma.masked_where(detectionSpindlesR==0, detectionSpindlesR==1)*(-0.032)
SOMaskedR=np.ma.masked_where(detectionSOR==0, detectionSOR==1)*(0.032)

#Plot A
timeA=np.arange(initT/fs,lenT/fs,1/fs)
timeC=np.arange(23,28,1/fs)
timeD=np.arange(14,17,1/fs)
axA1.plot(timeA,Ee[initT:lenT],color='black',label='SHAM')
axA1.plot(timeA,spindlesMasked[initT:lenT],'.-',ms=3,markevery=50,color=colorsProp[2])
axA1.plot(timeA,SOMasked[initT:lenT],'.-',ms=3,markevery=50,color=colorsProp[1])
axA1.text(-0.2,0.95,'A',fontsize=10,fontweight=1000,verticalalignment='bottom',transform=axA1.transAxes)

axA4.hlines(y=-0.03, xmin=23, xmax=28, color='gray', linewidth=0.5)
axA4.hlines(y=0.03, xmin=23, xmax=28, color='gray', linewidth=0.5)
axA4.vlines(x=23,ymin=-0.03,ymax=0.03,color='gray', linewidth=0.5)
axA4.vlines(x=28,ymin=-0.03,ymax=0.03,color='gray', linewidth=0.5)
axA4.plot(timeA,Phieph1[initT:lenT],color=colorsa[2],label='STIM-CL')
axA4.vlines(x=np.argwhere(markerph[initT:lenT]==1)/fs+initT/fs,ymin=0,ymax=offset+0.05,linewidth=1,color='gray')
axA4.plot(timeA,spindlesMaskedP[initT:lenT],'.-',ms=3,markevery=50,label=r'Spindles',color=colorsProp[2])
axA4.plot(timeA,SOMaskedP[initT:lenT],'.-',label='SO',ms=3,markevery=50,color=colorsProp[1])

axA3.vlines(x=np.argwhere(marker[initT:lenT]==1)/fs+initT/fs,ymin=0,ymax=offset+0.05,linewidth=1,color='gray')
axA3.plot(timeA,Phie1[initT:lenT],color=colorsa[0],label='STIM-P')
axA3.plot(timeA,spindlesMaskedC[initT:lenT],'.-',ms=3,markevery=50,color=colorsProp[2])
axA3.plot(timeA,SOMaskedC[initT:lenT],'.-',ms=3,markevery=50,color=colorsProp[1])

axA2.vlines(x=np.argwhere(markers[initT:lenT]==1)/fs+initT/fs,ymin=0,ymax=offset+0.05,linewidth=1,color='gray')
axA2.plot(timeA,Phiep1[initT:lenT],color=colorsa[1],label='STIM-R')
axA2.plot(timeA,spindlesMaskedR[initT:lenT],'.-',ms=3,markevery=50,color=colorsProp[2])
axA2.plot(timeA,SOMaskedR[initT:lenT],'.-',ms=3,markevery=50,color=colorsProp[1])



axA1.set_xlim([initT/fs, (lenT+1)/fs])
axA2.set_xlim([initT/fs, (lenT+1)/fs])
axA3.set_xlim([initT/fs, (lenT+1)/fs])
axA4.set_xlim([initT/fs, (lenT+1)/fs])
axA4.set_xlabel('Time (s)',fontsize=8,labelpad=0.3)
axA1.set_ylabel('$b(t)\ (s^{-1})$',fontsize=8,labelpad=0.05)
axA2.set_ylabel('$x(t)\ (s^{-1})$',fontsize=8,labelpad=0.05)
axA3.set_ylabel('$x(t)\ (s^{-1})$',fontsize=8,labelpad=0.05)
axA4.set_ylabel('$x(t)\ (s^{-1})$',fontsize=8,labelpad=0.05)
axA4.set_xticks(np.arange(initT/fs,lenT/fs+1,3))
# axA4.set_xticklabels(np.arange(initT/fs,lenT/fs,2),fontsize=8)
axA4.set_xticklabels(['15','18','21','24','27','30'],fontsize=8)
axA1.set_yticks([-0.04,0,0.04])
axA1.set_yticklabels(['-0.04','0.00','0.04'],fontsize=8)
axA2.set_yticks([-0.04,0,0.04])
axA2.set_yticklabels(['-0.04','0.00','0.04'],fontsize=8)
axA3.set_yticks([-0.04,0,0.04])
axA3.set_yticklabels(['-0.04','0.00','0.04'],fontsize=8)
axA4.set_yticks([-0.04,0,0.04])
axA4.set_yticklabels(['-0.04','0.00','0.04'],fontsize=8)
axA1.set_ylim([-0.04, 0.04])
axA2.set_ylim([-0.04, 0.04])
axA3.set_ylim([-0.04, 0.04])
axA4.set_ylim([-0.04, 0.04])

# legA1=axA1.legend(ncol=3,loc='upper left',fontsize=8,bbox_to_anchor=(-0.15,1.09,1,0.1))
axA1.tick_params(axis='y', which='both', labelsize=8)
axA2.tick_params(axis='y', which='both', labelsize=8)
axA3.tick_params(axis='y', which='both', labelsize=8)
axA4.tick_params(axis='y', which='both', labelsize=8)

axA1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
axA2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
axA3.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
axA4.tick_params(axis='x', which='both', labelsize=8)


#Plot C
axC1.plot(timeC,Phiephz[2300:2800],'b',label='$x(t)$')
axC1.plot(timeC,spindlesMaskedP[2300:2800]-3,'.-',ms=3,markevery=50,color=colorsProp[2],label=r'SP')
axC1.plot(timeC,SOMaskedP[2300:2800]+4,'.-',ms=3,markevery=50,color=colorsProp[1],label=r'SO')
axC1.vlines(x=np.argwhere(markerph[2300:2800]==1)/fs+2300/fs,ymin=-3,ymax=3,linewidth=0.8,color='gray')
axC1.text(-0.2,1.05,'C',fontsize=10,fontweight=1000,verticalalignment='bottom',transform=axC1.transAxes)
axC1.set_ylabel('z-score')
axC1.set_ylim([-4,5])
axC1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
axC1.tick_params(axis='y', which='both', labelsize=8)
axC1.legend(fontsize=8,ncol=3,loc='upper left',bbox_to_anchor=(-0.01,1.42,0.1,0.1),columnspacing=5.5,handletextpad=1.0)
#Filtered Bands
axC2.plot(timeC,filtered_sp[2300:2800],color=colorsProp[2],label='SP band')
axC2.plot(timeC,filtered_d[2300:2800],color=colorsProp[1],label='SO band')
axC2.plot(timeC,filtered_85[2300:2800]*60,'k',label='0.85 Hz')
axC2.vlines(x=np.argwhere(markerph[2300:2800]==1)/fs+2300/fs,ymin=-3,ymax=3,linewidth=0.8,color='gray')
# axC2.text(-0.25,0.88,'C2',fontsize=10,fontweight=1000,verticalalignment='bottom',transform=axC2.transAxes)
axC2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
axC2.tick_params(axis='y', which='both', labelsize=8)
axC2.set_ylabel('z-score')
axC2.legend(fontsize=8,ncol=3,loc='upper left',bbox_to_anchor=(-0.01,1.32,0.1,0.1),columnspacing=2.25,handletextpad=1.0)
axC2.set_yticks([-2,0,2])
axC2.set_yticklabels(['-2','0','2'])
#Phase
phase_offline=np.angle(signal.hilbert(filtered_85)*np.exp(-1j*np.pi/2))+np.pi

axC3.plot(timeC,phase[2300:2800]*0.3,'k',label='Online phase')
axC3.plot(timeC,phase_offline[2300:2800]*0.3,color=plt.cm.Greys(0.6),label='Offline phase')
axC3.vlines(x=np.argwhere(markerph[2300:2800]==1)/fs+2300/fs,ymin=-0.1,ymax=0.6*np.pi,linewidth=0.8,color='gray')
axC3.legend(fontsize=8,ncol=2,loc='upper left',bbox_to_anchor=(-0.01,1.42,0.1,0.1),columnspacing=7.2,handletextpad=1.0)
# axC3.text(-0.25,0.88,'C3',fontsize=10,fontweight=1000,verticalalignment='bottom',transform=axC3.transAxes)
axC3.set_ylabel('radians')
axC3.set_yticks([0,0.6*np.pi])
axC3.set_ylim([-0.1,0.65*np.pi])
axC3.set_xticks([23,24,25,26,27,28])
axC3.set_xticklabels([23,24,25,26,27,28])
axC3.set_yticklabels(['0',r'$2\pi$'],fontsize=8)
axC3.set_xlabel('Time (s)',fontsize=8,labelpad=0.2)
axC3.tick_params(axis='both', which='both', labelsize=8)


#Plot B
normSW=colors.PowerNorm(gamma=0.5,vmin=0,vmax=1)

imB1=axB1.imshow((np.abs(ScalogramB[33:216,initT:lenT])),extent=(-0.5,lenT+0.5,32.5,215.5),aspect='auto',norm=normSW,cmap=cm.CMRmap)
axB1.set_ylabel('Freq. (Hz)',fontsize=8)
axB1.set_yticks(indexesWavelet)
axB1.set_yticklabels(labels_appFreq,fontsize=8)
axB1.set_xlim([initT, lenT+1])
axB1.text(-0.28,0.95,'B',fontsize=10,fontweight=1000,verticalalignment='bottom',transform=axB1.transAxes)
axB1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

imB2=axB2.imshow((np.abs(ScalogramRandom[33:216,initT:lenT])),extent=(-0.5,lenT+0.5,32.5,215.5),aspect='auto',norm=normSW,cmap=cm.CMRmap)
axB2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
axB2.set_ylabel('Freq. (Hz)',fontsize=8)
axB2.set_yticks(indexesWavelet)
axB2.set_yticklabels(labels_appFreq,fontsize=8)
axB2.set_xlim([initT, lenT+1])


imB3=axB3.imshow((np.abs(Scalogram[33:216,initT:lenT])),extent=(-0.5,lenT+0.5,32.5,215.5),aspect='auto',norm=normSW,cmap=cm.CMRmap)
axB3.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
axB3.set_ylabel('Freq. (Hz)',fontsize=8)
axB3.set_yticks(indexesWavelet)
axB3.set_yticklabels(labels_appFreq,fontsize=8)
axB3.set_xlim([initT, lenT+1])


imB4=axB4.imshow((np.abs(ScalogramPhase[33:216,initT:lenT])),extent=(-0.5,lenT+0.5,32.5,215.5),aspect='auto',norm=normSW,cmap=cm.CMRmap)
axB4.set_xticks(np.arange(initT,lenT+fs,3*fs))
axB4.set_xticklabels(['15','18','21','24','27','30'],fontsize=8)
axB4.set_ylabel('Freq. (Hz)',fontsize=8)
axB4.set_yticks(indexesWavelet)
axB4.set_yticklabels(labels_appFreq,fontsize=8)
axB4.set_xlabel('Time (s)',fontsize=8,labelpad=0.2)
axB4.set_xlim([initT, lenT+1])


clb=fig1.colorbar(imB4,cax=axB5,ax=axB5,orientation='vertical',fraction=80,shrink=1.0)
clb.set_label('(a. u.)',fontsize=8,labelpad=-38)
axB5.yaxis.set_ticks_position('left')
axB5.tick_params(axis='both', which='both', labelsize=8)

##
#Zoom events

# axD.plot(timeD,Phieph1[1400:1700],'b',label='$x(t)$')
# axD.plot(timeD,spindlesMaskedP[1400:1700]+0.02,'.-',ms=3,markevery=50,color=colorsProp[4],label=r'SP')
# axD.plot(timeD,SOMaskedP[1400:1700]-0.002,'.-',ms=3,markevery=50,color=colorsProp[5],label=r'SO')
# axD.set_ylabel('x(t) $s^{-1}$',fontsize=8,labelpad=0.12)
# axD.set_xlabel('Time (s)',fontsize=8,labelpad=0.2)
# axD.tick_params(axis='both', which='both', labelsize=8)
# axD.legend(fontsize=8,ncol=3,loc='upper left',bbox_to_anchor=(-0.01,1.43,0.1,0.1),columnspacing=4.0,handletextpad=1.0)
# axD.text(-0.2,1.00,'D',fontsize=10,fontweight=1000,verticalalignment='bottom',transform=axD.transAxes)

##
#Sum in time
ScalogramB=np.load('../collectedData/ScalogramB.npz')['ScalogramB']
ScalogramP=np.load('../collectedData/ScalogramP.npz')['ScalogramP']
ScalogramR=np.load('../collectedData/ScalogramR.npz')['ScalogramR']
ScalogramCL=np.load('../collectedData/ScalogramCL.npz')['ScalogramCL']

axE.loglog(freqsWavelet,ScalogramB/91000,'k',label='SHAM')
axE.loglog(freqsWavelet,ScalogramR/91000,color=colorsa[1],linewidth=1.5,label='STIM-R')
axE.loglog(freqsWavelet,ScalogramP/91000,color=colorsa[0],linewidth=1.5,label='STIM-P')
axE.loglog(freqsWavelet,ScalogramCL/91000,color=colorsa[2],linewidth=1.5,label='STIM-CL 45')

axE.fill_between(freqsWavelet[166:214],ScalogramB[166:214]/91000,ScalogramCL[166:214]/91000,color='#666666')
axE.fill_between(freqsWavelet[166:214],np.ones((48,))*1e-3,ScalogramB[166:214]/91000,color='#aaaaaa')
axE.fill_between(freqsWavelet[33:65],np.ones((32,))*1e-3,ScalogramB[33:65]/91000,color='#aaaaaa')
axE.fill_between(freqsWavelet[33:65],ScalogramB[33:65]/91000,ScalogramCL[33:65]/91000,color='#777777')

axE.text(0.7,0.3,'SO',fontsize=9)
axE.text(11,0.105,'SP',fontsize=9)
axE.set_xlim([0.5,20.5])

axE.set_xlabel('Frequency (Hz)',fontsize=8,labelpad=-0.1)
axE.set_ylabel('a. u.',fontsize=8,labelpad=0)
axE.set_ylim([1e-2,0.8])
axE.text(-0.2,1.05,'D',fontsize=10,fontweight=1000,verticalalignment='bottom',transform=axE.transAxes)
axE.tick_params(axis='both', which='both', labelsize=8)

#Analytical spectrum
print("Calculate spectrum with phi_n=1")
matfile=h5py.File('/home/felipe/Dropbox/UTFSM/NeuralField/analyze/collectedData/NF-interpolation/interpolation_step_66-phin1.mat','r')
solutions_e=np.array(matfile['soulutions_e'])
phin=1
index_phin=np.argwhere(solutions_e[:,0]>=phin)[0][0]
phi_n=solutions_e[index_phin,0]
phi_e=solutions_e[index_phin,1]
phi_r=solutions_e[index_phin,2]
phi_s=solutions_e[index_phin,3]
#Model Parameters
alpha=45
beta=186
gamma=116
r=0.086
k0=1
t0=0.085
Lx=0.5
Ly=0.5
Qmax=340
sigma_p=0.0038
nus=np.array([5.541, -5.652, 1.53, 0.286, 1.12, 2.67, -1.73, 9.22])*1e-3
omega=np.arange(0,2*np.pi*50+0.001,np.pi/100)
#Calcualte rho
rho_e=phi_e/sigma_p*(1-(phi_e/Qmax));
rho_r=phi_r/sigma_p*(1-(phi_r/Qmax));
rho_s=phi_s/sigma_p*(1-(phi_s/Qmax));
Gee=nus[0]*rho_e
Gei=nus[1]*rho_e
Ges=nus[2]*rho_e
Gre=nus[3]*rho_r
Grs=nus[4]*rho_r
Gse=nus[5]*rho_s
Gsr=nus[6]*rho_s
Gsn=nus[7]*rho_s
Gesre=Ges*Gsr*Gre
Gese=Ges*Gse
Gsrs=Gsr*Grs
rho=np.array([[rho_e],[rho_r],[rho_s]])
XYZ=np.array([Gee/(1-Gei), (Ges*Gse+Ges*Gsr*Gre)/((1-Gsr*Grs)*(1-Gei)),-Gsr*Grs*(alpha*beta)/(alpha+beta)**2])
Power=steadyFunctions.SpatialSpectrum(omega,Lx,Ly,alpha,beta,gamma,r,rho,nus,t0,k0)
#Correct power by frequency
correctedPower=np.zeros((np.shape(Power)[0]-1,))
correctedPower=Power[1::]*omega[1::]/(2*np.pi)
axE.loglog(omega[1::]/(2*np.pi),correctedPower*np.max(ScalogramB/91000)/np.max(correctedPower),'--r',label='Analytical')
axE.legend(fontsize=8,loc='upper left',ncol=3,bbox_to_anchor=(-0.005,0.91,1,0.1),columnspacing=2.7,handletextpad=1.0)


axF.hist([np.hstack([durationsSPB,durationsSPP,durationsSPR]),np.hstack([durationsSOB,durationsSOP,durationsSOR])],bins=np.arange(0.05,2,0.1),color=[colorsProp[2],colorsProp[1]],density=False)
axG.hist([np.hstack([p2pSPB,p2pSPP,p2pSPR]),np.hstack([p2pSOB,p2pSOP,p2pSOR])],bins=np.arange(0.0025,0.06,0.005),color=[colorsProp[2],colorsProp[1]],density=False)
axF.legend(['SP','SO'],loc='upper left',fontsize=8,handletextpad=0.2)
axG.legend(['SP','SO'],loc='upper left',fontsize=8,handletextpad=0.2)
axF.set_xlabel('Duration (s)',fontsize=8)
axF.set_ylabel('\# of events',fontsize=8)
axG.set_xlabel('Peak-to-peak Amplitude ($s^{-1}$)',fontsize=8)
axG.set_ylabel('\# of events',fontsize=8)
axF.tick_params(axis='both', which='both', labelsize=8)
axG.tick_params(axis='both', which='both', labelsize=8)
axF.text(-0.19,1.05,'E',fontsize=10,fontweight=1000,verticalalignment='bottom',transform=axF.transAxes)
axG.text(-0.19,1.05,'F',fontsize=10,fontweight=1000,verticalalignment='bottom',transform=axG.transAxes)



fig1.savefig('./output/Fig2.eps',dpi=300,bbox_inches='tight',bbox_extra_artists=[axG,axB5])