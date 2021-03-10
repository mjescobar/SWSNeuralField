# -*- coding: utf-8 -*-
"""
Created on Tue Oct 04 09:58:47 2016

@author: porio
"""


from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import Wavelets

dt=0.02
time=np.arange(0,100,dt)
#x1+=np.random.normal(size=len(x1))*0.4

x=-2 + 6*np.sin(2*np.pi*5*time) +5*np.sin(2*np.pi*2*time)
#x+=np.random.normal(size=len(x))*0.8

x2=1*np.cos(2*np.pi*(time-39)**1.5) * (time-40)*np.exp(-(time-39)/10)
x2[(time<40)+(time>68)]=0
x+=x2
   
plt.figure(4)
plt.clf()
plt.plot(time,x,alpha=1)
plt.xlabel('Time (s)')

#%%
fftx=np.fft.fft(x-np.mean(x))

ffreq=np.fft.fftfreq(len(x),dt)

limplot=range(1001)

plt.figure(3)
plt.clf()
plt.plot(ffreq[limplot],np.abs(fftx[limplot])**2)
plt.xlim((0,10))
plt.yscale('log');plt.ylim((0.1,0.5e8))

plt.xlabel('Freq (Hz)')

#%%
plt.figure(2)
plt.clf()
spectrum,freqs,t,_=plt.specgram(x-np.mean(x),Fs=50,NFFT=128,noverlap=64)
#spectrum,freqs,t,_=plt.specgram(x1-np.mean(x1),Fs=50)
plt.xlabel('Time (s)')
plt.ylabel('Frequency (Hz)')


#%%
# freqs = np.arange(0.1,15,0.1)  # Hz
freqs = np.logspace(-0.3,1.47,150)
# Periods = 1/(freqs*dt)    #Desired periods in sample untis. dt has to be in seconds
# dScales=Periods/Wavelets.Morlet.fourierwl  #desired Scales

wavel1=Wavelets.Morlet(x,freqs=freqs,dt=dt,omega0=10)

scales=wavel1.getscales()

cwt1=wavel1.getdata()
cwt1=cwt1/(np.sqrt(scales[:,None]))
#pwr1=wavel1.getpower()#/(np.sqrt(scales[:,None]))
#pwr1=np.sqrt(np.real(cwt1*np.conjugate(cwt1)))
pwr1=wavel1.getnormpower()

fmin=min(freqs)  #max/min Frequencies in Hz
fmax=max(freqs)
#%%
plt.figure(1)
plt.clf()

ax1=plt.subplot2grid((1,5),(0,0),colspan=4)
#
# plt.imshow(pwr1,cmap='jet',extent=(min(t),max(t),fmin,fmax),origin='lower', #Frequency decreases from top to bottom
#            interpolation='none',aspect='auto')

plt.imshow(pwr1,cmap='jet',extent=(min(t),max(t),np.log10(fmin),np.log10(fmax)),origin='lower', #Frequency decreases from top to bottom
           interpolation='none',aspect='auto')

# ax1.set_yscale('log')
ax1.set_ylabel('Frequency (Hz)')
##
#plt.imshow(np.absolute(pwr1),cmap='jet',
#           extent=(min(t),max(t),Tmin,Tmax),origin='lower',  #Period increases from bottom to top
#           interpolation='none',aspect='auto')
#ax1.set_yscale('log')
#ax1.set_ylabel('Period (s)')

#ax4=plt.subplot2grid((5,5),(1,4),rowspan=2)
#plt.plot(np.sum(pwr,axis=1),fourierwl*dt)
#ax4.set_yscale('log')
#ax4.set_ylim((Tmin,Tmax))

ax5=plt.subplot2grid((1,5),(0,4),rowspan=2)
plt.plot(np.sum(pwr1,axis=1),freqs)  #plotted against frequency, 1/period
#plt.plot(np.sum(pwr1,axis=1),fourierwl*dt)  #plotted against frequency, 1/period
ax5.set_yscale('log')
ax5.set_ylim((fmin,fmax))
#ax5.set_ylim((Tmin,Tmax))

plt.tight_layout()


