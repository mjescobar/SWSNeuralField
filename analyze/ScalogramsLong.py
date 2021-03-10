#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 20:03:47 2021

@author: felipe
"""

import sys, os
sys.path.append(os.path.abspath(os.path.join('.', 'spectrum')))
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
import ERP
import Wavelets
import re

fileScalogramP='/media/felipe/TOSHIBA2T/Frequency-Long/SWS-seed-1-decreaseRamp-F0.85-A3.46e+01-D1.00e-01-sdA-OFF-sdF-OFF-pF0.02-pR2.00-Phi.txt'
fileScalogramR='/media/felipe/TOSHIBA2T/Frequency-Long/SWS-seed-1-decreaseRamp-F0.85-A3.46e+01-D1.00e-01-sdA-OFF-sdF-Poisson-pF0.02-pR2.00-Phi.txt'
fileScalogramCL='/media/felipe/TOSHIBA2T/Phase-Long/SWSlong-seed-1-phase-0.79-decreaseRamp-F0.85-A3.46e+01-D1.00e-01-sdA-OFF-sdF-OFF-pF0.02-pR2.00-Phi.txt'
fileScalogramBase='/media/felipe/TOSHIBA2T/Baseline-Long/SWSlong-seed-1-baseline-Phi.txt'

PhieS1,PhieSz,PhinSz,BaselineeSz,BaselineSn,StimS,markerS,timeS=Data.loadAllData(fileScalogramP, fileScalogramBase)
PhieR1,PhieRz,PhinRz,BaselineeRz,BaselineRn,StimR,markerR,timeR=Data.loadAllData(fileScalogramR, fileScalogramBase)
PhiephS1,PhiephSz,PhinphSz,BaselineephSz,BaselinephSn,StimphS,markerphS,timephS=Data.loadAllData(fileScalogramCL, fileScalogramBase)
print('ScalogramP')
freqsWavelet,scales,Scalogram=Data.waveletSingle(PhieSz,correctF=False)
np.savez('ScalogramP.npz',ScalogramP=np.sum(Scalogram,axis=1))
del freqsWavelet,scales,Scalogram
print('ScalogramR')
freqsWavelet,scales,ScalogramRandom=Data.waveletSingle(PhieRz,correctF=False)
np.savez('ScalogramR.npz',ScalogramR=np.sum(ScalogramRandom,axis=1))
del freqsWavelet,scales,ScalogramRandom
print('ScalogramB')
freqB,scalesB,ScalogramB=Data.waveletSingle(BaselineeSz,correctF=False)
np.savez('ScalogramB.npz',ScalogramB=np.sum(ScalogramB,axis=1))
del freqB,scalesB,ScalogramB
print('ScalogramCL')
freqsPhase,scalesPhase,ScalogramPhase=Data.waveletSingle(PhiephSz,correctF=False)
np.savez('ScalogramCL.npz',ScalogramCL=np.sum(ScalogramPhase,axis=1))
del freqsPhase,scalesPhase,ScalogramPhase
#%%



