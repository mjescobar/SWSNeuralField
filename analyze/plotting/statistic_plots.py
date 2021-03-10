#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 23:24:04 2020

@author: felipe
"""
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.axes as axes

def polar_statistics(angles,axis=0,nan_policy='omit'):
    #return the mean and standard deviation of circular data
    #np.circmean returns same value 
    #np.circstd returns different value
    sin_angles=np.sin(angles)
    cos_angles=np.cos(angles)
    meanpolar=np.arctan(np.mean(sin_angles,axis=axis)/np.mean(cos_angles,axis=axis))
    stdpolar=np.arctan(np.std(sin_angles,axis=axis)/np.std(cos_angles,axis=axis))
    
    return meanpolar,stdpolar
def circ_hist(a,ax,bins=10,bottom=1,mean_axis=0,density=False,cmap=plt.cm.Blues,plotmeanstd=True):
    #Circular histogram with mean value and standard deviation
    #a:angles data
    #ax: matplotlib.pyplot.axes with projection='polar'
    #bins: number of bins in the range [0,360) or [0,2pi)
    #mean_axis: mean calculated along this axis
    #density: pmf type histogram
    #cmap: colorma
    #plotmeanstd: plot or not the mean and standard deviation markers
    #return the matplotlib.pyplot.axes
    ###Calcualte the histogram
    hist,bin_edges=np.histogram(a,bins=bins,density=density)
    width = (2*np.pi) / (1.2*bins)
    ax = ax
    ###Calculate the statistics
    mean_theta=stats.circmean(a,axis=mean_axis,nan_policy='omit')
    std_theta=stats.circstd(a,axis=mean_axis,nan_policy='omit')
    ###Create circular bars plot from the histogram data
    bars = ax.bar(bin_edges[0:-1], hist, width=width, bottom=bottom)
    ###Plot the statisitics markers
    if plotmeanstd:
        ax.bar(mean_theta, np.max(hist), width=width/3, bottom=bottom, facecolor='red')
        ax.bar(mean_theta+std_theta, np.max(hist), width=width/4, bottom=bottom, facecolor='black')
        ax.bar(mean_theta-std_theta, np.max(hist), width=width/4, bottom=bottom, facecolor='black')
    ##Change the color of the bars
    # Use custom colors and opacity
    for r, bar in zip(hist, bars):
        bar.set_facecolor(cmap(r / 10.))
    return ax


def bivariable_std(datax,datay,ax,meanx_zero=0,meany_zero=0,stdx_zero=1,stdy_zero=1,mean_axis=0,cmap=plt.cm.Blues,marker='o',
                   flagstdx=False,flagstdy=False,meanzero_type='dot',relative_diference=True,labels=[''],colormap='None'):
    #datax and datay: 1-D array
    #meanx_zero and meany_zero: means of baseline data
    #stdx_zero and stdy_zero: standard deviation of baseline data
    #cmap: colormap
    #marker: type of marker
    #flagstdx: plot the standard deviation in the x-axis
    #flagstdy: plot the stadard deviation in the y-axis
    #relative_diference: plot in percentages of realive error with respect the means of the baseline data
    #return ax, a matplotlib.pyplot.axes instance
    ###Statisitic of the data
    flag_error=0
    if datay=='None':
        meanx=np.nanmean(datax,axis=mean_axis)
        stdx=np.nanstd(datax,axis=mean_axis)
        stdy=np.zeros_like(stdx)
        meany=np.arange(len(meanx))
        flagstdy=False
    elif datax=='None':
        meany=np.nanmean(datay,axis=mean_axis)
        stdy=np.nanstd(datay,axis=mean_axis)
        stdx=np.zeros_like(stdy)
        meanx=np.arange(len(meany))
        flagstdx=False
    else:
        meanx=np.nanmean(datax,axis=mean_axis)
        stdx=np.nanstd(datax,axis=mean_axis)
        meany=np.nanmean(datay,axis=mean_axis)
        stdy=np.nanstd(datay,axis=mean_axis)
    
    if type(meanx)!=np.float64:
        data_length=len(meanx)
        if len(meanx)!=len(meany):
            flag_error=2
    else:
        data_length=1
    ###If relative difference is required
    if relative_diference:
        ##If means is zero, the increase or decrease is directly a percentage/
        if meanx_zero==0:
            meanx=meanx*100 
            stdx=stdx*100
        else:
            meanx=(meanx-meanx_zero)/(meanx_zero)*100
            stdx=(stdx)/(meanx_zero)*100
        
        if meany_zero==0:
            meany=meany*100
            stdx=stdx*100
        else:
            meany=(meany-meany_zero)/(meany_zero)*100
            stdy=(stdy)/(meany_zero)*100
        #Change baseline values to zero
        if meanx_zero==0:
            stdx_zero=stdx_zero*100
        else:
            stdx_zero=(stdx_zero)/(meanx_zero)*100
        
        if meany_zero==0:
            stdy_zero=stdy_zero*100
        else:
            stdy_zero=(stdy_zero)/(meany_zero)*100
        meanx_zero=0
        meany_zero=0
       
    
    if flag_error==0:
        if meanzero_type=='dot':
            ax.plot(meanx_zero,meany_zero,'k',marker='o',markersize=10,label='SHAM')
            if flagstdx:
                ax.plot([stdx_zero+meanx_zero,stdx_zero+meanx_zero],[np.min([meany_zero-stdy_zero,np.min(meany-stdy)]),np.max([meany_zero+stdy_zero,np.max(meany+stdy)])],'+k')
                ax.plot([meanx_zero-stdx_zero,meanx_zero-stdx_zero],[np.min([meany_zero-stdy_zero,np.min(meany-stdy)]),np.max([meany_zero+stdy_zero,np.max(meany+stdy)])],'+k')
            if flagstdy:
                ax.plot([np.min([meanx_zero-stdx_zero,np.min(meanx-stdx)]),np.max([meanx_zero+stdx_zero,np.max(meanx+stdx)])],[stdy_zero+meany_zero,stdy_zero+meany_zero],'+k')
                ax.plot([np.min([meanx_zero-stdx_zero,np.min(meanx-stdx)]),np.max([meanx_zero+stdx_zero,np.max(meanx+stdx)])],[meany_zero-stdy_zero,meany_zero-stdy_zero],'+k')
        elif meanzero_type=='line':
            ax.plot([meanx_zero,meanx_zero],[np.min([meany_zero-stdy_zero,np.min(meany-stdy)]),np.max([meany_zero+stdy_zero,np.max(meany+stdy)])],'k',linewidth=0.5,label='SHAM')
            ax.plot([np.min([meanx_zero-stdx_zero,np.min(meanx-stdx)]),np.max([meanx_zero+stdx_zero,np.max(meanx+stdx)])],[meany_zero,meany_zero],'k',linewidth=0.5)
            if flagstdx:
                ax.plot([stdx_zero+meanx_zero,stdx_zero+meanx_zero],[np.min([meany_zero-stdy_zero,np.min(meany-stdy)]),np.max([meany_zero+stdy_zero,np.max(meany+stdy)])],'k',linestyle=(0,(2,5)),linewidth=0.5)
                ax.plot([meanx_zero-stdx_zero,meanx_zero-stdx_zero],[np.min([meany_zero-stdy_zero,np.min(meany-stdy)]),np.max([meany_zero+stdy_zero,np.max(meany+stdy)])],'k',linestyle=(0,(2,5)),linewidth=0.5)
            if flagstdy:
                ax.plot([np.min([meanx_zero-stdx_zero,np.min(meanx-stdx)]),np.max([meanx_zero+stdx_zero,np.max(meanx+stdx)])],[stdy_zero+meany_zero,stdy_zero+meany_zero],'k',linestyle=(0,(2,5)),linewidth=0.5)
                ax.plot([np.min([meanx_zero-stdx_zero,np.min(meanx-stdx)]),np.max([meanx_zero+stdx_zero,np.max(meanx+stdx)])],[meany_zero-stdy_zero,meany_zero-stdy_zero],'k',linestyle=(0,(2,5)),linewidth=0.5)
        elif meanzero_type=='line_inf':
            ax.plot([meanx_zero,meanx_zero],[-500,1000],'k',linewidth=0.5,label='SHAM')
            ax.plot([-500,1000],[meany_zero,meany_zero],'k',linewidth=0.5)
            if flagstdx:
                ax.plot([stdx_zero+meanx_zero,stdx_zero+meanx_zero],[-500,1000],'k',linestyle=(0,(2,5)),linewidth=0.5)
                ax.plot([meanx_zero-stdx_zero,meanx_zero-stdx_zero],[-500,1000],'k',linestyle=(0,(2,5)),linewidth=0.5)
            if flagstdy:
                ax.plot([-500,1000],[stdy_zero+meany_zero,stdy_zero+meany_zero],'k',linestyle=(0,(2,5)),linewidth=0.5)
                ax.plot([-500,1000],[meany_zero-stdy_zero,meany_zero-stdy_zero],'k',linestyle=(0,(2,5)),linewidth=0.5)

        elif meanzero_type=='None':
            print('Not plot of the mean')

        for x,y,sx,sy,n in zip(meanx,meany,stdx,stdy,range(data_length)):
            #Multiple values (Matrix data)
            if len(labels)>1:
                if colormap=='direct':
                    ax.plot(x,y,marker=marker,color=cmap(n),markersize=6,label=labels[n])
                    if flagstdx:
                        ax.plot([x-sx,x+sx],[y,y],':',color=cmap(n),linewidth=0.5)
                    if flagstdy:
                        ax.plot([x,x],[y-sy,y+sy],':',color=cmap(n),linewidth=0.5)
                else:
                    ax.plot(x,y,marker=marker,color=cmap(1.0-n*.12),markersize=6,label=labels[n]) 
                    if flagstdx:
                        ax.plot([x-sx,x+sx],[y,y],':',color=cmap(0.8),linewidth=0.5)
                    if flagstdy:
                        ax.plot([x,x],[y-sy,y+sy],':',color=cmap(0.6),linewidth=0.5)
            else:
                #Single value (array data; shape=(N,1))
                ax.plot(x,y,marker=marker,color=cmap,markersize=6,label=labels[0])
                if flagstdx:
                    ax.plot([x-sx,x+sx],[y,y],':',color=cmap,linewidth=0.5)
                if flagstdy:
                    ax.plot([x,x],[y-sy,y+sy],':',color=cmap,linewidth=0.5)
    else:
        print("Error in the length of data,len(datax) should be equal to len(datay)")
    return ax


def graph_l(dataxSO,dataySO,dataxSP,dataySP,dataPSO,dataPCP,E,x_axis,ax,mean_axis=0,cmap=plt.cm.Blues,marker='o',colormap='None'):
    meanxSO=np.nanmean(dataxSO,axis=mean_axis)
    meanySO=np.nanmean(dataySO,axis=mean_axis)
    meanxSP=np.nanmean(dataxSP,axis=mean_axis)
    meanySP=np.nanmean(dataySP,axis=mean_axis)
    meanPSO=np.nanmean(dataPSO,axis=mean_axis)
    meanPCP=np.nanmean(dataPCP,axis=mean_axis)
    if type(meanxSO)!=np.float64:
        data_length=len(meanxSO)
        if len(meanxSO)!=len(meanySO):
            flag_error=2
    else:
        data_length=1        

    l_all=(meanPCP-meanPSO)*(meanxSO*meanySO+meanxSP*meanySP)-E
    for x,l,n in  zip(x_axis,l_all,range(data_length)):
        if colormap=='direct':
            ax.plot(x,l,marker=marker,color=cmap(n))
        else:
            ax.plot(x,l,marker=marker,color=cmap(1.0-n*.12))
    return ax
    
def univariable_std(datax,datay,ax,meanx_zero=0,meany_zero=0,stdx_zero=1,stdy_zero=1,mean_axis=0,cmap=plt.cm.Blues,marker='o',
                   flagstdx=False,flagstdy=False,meanzero_type='dot',relative_diference=True,labels=[''],colormap='None',bt=3):
    #datax and datay: 1-D array
    #meanx_zero and meany_zero: means of baseline data
    #stdx_zero and stdy_zero: standard deviation of baseline data
    #cmap: colormap
    #marker: type of marker
    #flagstdx: plot the standard deviation in the x-axis
    #flagstdy: plot the stadard deviation in the y-axis
    #relative_diference: plot in percentages of realive error with respect the means of the baseline data
    #return ax, a matplotlib.pyplot.axes instance
    ###Statisitic of the data
    flag_error=0
    # ax1 = ax.twinx()
    if datay=='None':
        meanx=np.nanmean(datax,axis=mean_axis)
        stdx=np.nanstd(datax,axis=mean_axis)
        stdy=np.zeros_like(stdx)
        meany=np.arange(len(meanx))
        flagstdy=False
    elif datax=='None':
        meany=np.nanmean(datay,axis=mean_axis)
        stdy=np.nanstd(datay,axis=mean_axis)
        stdx=np.zeros_like(stdy)
        meanx=np.arange(len(meany))
        flagstdx=False
    else:
        meanx=np.nanmean(datax,axis=mean_axis)
        stdx=np.nanstd(datax,axis=mean_axis)
        meany=np.nanmean(datay,axis=mean_axis)
        stdy=np.nanstd(datay,axis=mean_axis)
    
    if type(meanx)!=np.float64:
        data_length=len(meanx)
        if len(meanx)!=len(meany):
            flag_error=2
    else:
        data_length=1
    ###If relative difference is required
    if relative_diference:
        ##If means is zero, the increase or decrease is directly a percentage/
        if meanx_zero==0:
            meanx=meanx*100 
            stdx=stdx*100
        else:
            meanx=(meanx-meanx_zero)/(meanx_zero)*100
            stdx=(stdx)/(meanx_zero)*100
        
        if meany_zero==0:
            meany=meany*100
            stdx=stdx*100
        else:
            meany=(meany-meany_zero)/(meany_zero)*100
            stdy=(stdy)/(meany_zero)*100
        #Change baseline values to zero
        if meanx_zero==0:
            stdx_zero=stdx_zero*100
        else:
            stdx_zero=(stdx_zero)/(meanx_zero)*100
        
        if meany_zero==0:
            stdy_zero=stdy_zero*100
        else:
            stdy_zero=(stdy_zero)/(meany_zero)*100
        meanx_zero=0
        meany_zero=0
       
    
    if flag_error==0:
        if meanzero_type=='line':
            ax.plot([0,data_length],[meanx_zero,meanx_zero],'k',linewidth=2,label='SHAM')
            ax.plot([0,data_length],[meany_zero-bt,meany_zero-bt],color='gray',linewidth=2)
            if flagstdx:
                ax.plot([0,data_length],[stdx_zero+meanx_zero,stdx_zero+meanx_zero],'k',linestyle=(0,(2,5)))
                ax.plot([0,data_length],[meanx_zero-stdx_zero,meanx_zero-stdx_zero],'k',linestyle=(0,(2,5)))
            if flagstdy:
                ax.plot([0,data_length],[stdy_zero+meany_zero-bt,stdy_zero+meany_zero-bt],color='gray',linestyle=(0,(2,5)))
                ax.plot([0,data_length],[meany_zero-stdy_zero-bt,meany_zero-stdy_zero-bt],color='gray',linestyle=(0,(2,5)))
        elif meanzero_type=='None':
            print('Not plot of the mean')

        for x,y,sx,sy,n in zip(meanx,meany,stdx,stdy,range(data_length)):
            #Multiple values (Matrix data)
            if len(labels)>1:
                if colormap=='direct':
                    if flagstdx:
                        ax.bar(n,x,yerr=sx,color=cmap(n),tick_label=labels[n])
                    else:
                        ax.bar(n,x,color=cmap(n))
                    if flagstdy:
                        ax.bar(n,y-bt,yerr=sy,color=cmap(n),linewidth=2,bottom=bt,tick_label=labels[n])
                    else:
                        ax.bar(n,y-bt,color=cmap(n),bottom=bt,linewidth=2)
                else:
                    if flagstdx:
                        ax.bar(n,x,yerr=sx,color=cmap(1.0-n*.12),tick_label=labels[n]) 
                    else:
                        ax.bar(n,x,color=cmap(1.0-n*.12),tick_label=labels[n]) 
                    if flagstdy:
                        ax.bar(n,y-bt,yerr=sy,color=cmap(1.0-n*.12),bottom=bt,linewidth=2)
                    if flagstdy:
                        ax.bar(n,y-bt,color=cmap(1.0-n*.12),bottom=bt,linewidth=2) 
            else:
                #Single value (array data; shape=(N,1))
                ax.plot([0,data_length],[x,x],marker=marker,color=cmap,markersize=6,label=labels[0])
                ax.plot([0,data_length],[y-bt,y-bt],marker=marker,color=cmap,markersize=6,label=labels[0])
                if flagstdx:
                    ax.plot([0,data_length],[x-sx,x-sx],':',color=cmap)
                    ax.plot([0,data_length],[x+sx,x+sx],':',color=cmap)
                if flagstdy:
                    ax.plot([0,data_length],[y-sy-bt,y-sy-bt],':',color=cmap)
                    ax.plot([0,data_length],[y+sy-bt,y+sy-bt],':',color=cmap)
    else:
        print("Error in the length of data,len(datax) should be equal to len(datay)")
    return ax

def matrix_ttest(data,ax,test_type='ind',p_level=0.01,color='black',cmap=plt.cm.viridis,axis=0,labels=[''],annotate=False,equal_var=False):
    
    shape_matrix=np.shape(data)
    len_data=shape_matrix[axis]
    matrix_t=np.zeros((len_data,len_data))
    matrix_p=np.ones((len_data,len_data))
    for i in range(len_data):
        if axis==0:
            pret,prep=stats.shapiro(data[i,:])
            # print('i: ',pret,prep)
        else:
            pret,prep=stats.shapiro(data[:,i])
            # print('i: ',pret,prep)
        for j in range(len_data):
            if i>j:
                if axis==0:
                    pret,prepj=stats.shapiro(data[j,:])
                else:
                    pret,prepj=stats.shapiro(data[:,j])
                if prep<=0.05 or prepj<=0.05:
                    if axis==0:
                        matrix_t[i,j],matrix_p[i,j]=stats.wilcoxon(data[i,:],data[j,:],mode = 'approx')
                    else:
                        matrix_t[i,j],matrix_p[i,j]=stats.wilcoxon(data[:,i],data[:,j],mode = 'approx')
                else:    
                    if test_type=='ind':
                        if axis==0:
                            matrix_t[i,j],matrix_p[i,j]=stats.ttest_ind(data[i,:],data[j,:],equal_var=equal_var)
                        else:
                            matrix_t[i,j],matrix_p[i,j]=stats.ttest_ind(data[:,i],data[:,j],equal_var=equal_var)
                    elif test_type=='rel':
                        if axis==0:
                            matrix_t[i,j],matrix_p[i,j]=stats.ttest_rel(data[i,:],data[j,:])
                        else:
                            matrix_t[i,j],matrix_p[i,j]=stats.ttest_rel(data[:,i],data[:,j])
                # print('i:',i,'j:',j,'&','%.1f'%matrix_t[i,j],'& %.1e'%matrix_p[i,j])
    im=ax.imshow(matrix_t,aspect='auto',extent=(0,len_data,0,len_data),cmap=cmap)
    significant_pixels=np.argwhere(matrix_p<p_level)
    
    for pix in significant_pixels:
        i=pix[0]
        j=pix[1]
        ax.axhline(len_data-i,j/len_data,(j+1)/len_data,color=color)
        ax.axhline(len_data-i-0.995,j/len_data,(j+1)/len_data,color=color)
        ax.axvline(j+0.01,(len_data-i-1)/len_data,(len_data-i)/len_data,color=color)
        ax.axvline(j+0.995,(len_data-i-1)/len_data,(len_data-i)/len_data,color=color)
        ax.plot([j,j+1],[len_data-i-1,len_data-i],color=color)
        if annotate:
            ax.text(j+0.05,len_data-i-0.9,'%.1e'%matrix_p[i,j],fontsize=8)
    
    return im,ax

def vsmatrix_ttest(data,data1,ax,test_type='ind',p_level=0.01,color='black',cmap=plt.cm.viridis,axis=0,labels=[''],annotate=False,equal_var=False):
    shape_matrix=np.shape(data)
    len_data=shape_matrix[axis]
    matrix_t=np.zeros((len_data,len_data))
    matrix_p=np.ones((len_data,len_data))
    for i in range(len_data):
        for j in range(len_data):
                if test_type=='ind':
                    if axis==0:
                        matrix_t[i,j],matrix_p[i,j]=stats.ttest_ind(data[i,:],data1[j,:],equal_var=equal_var)
                    else:
                        matrix_t[i,j],matrix_p[i,j]=stats.ttest_ind(data[:,i],data1[:,j],equal_var=equal_var)
                elif test_type=='rel':
                    if axis==0:
                        matrix_t[i],matrix_p[i,j]=stats.ttest_rel(data[i,:],data1[j,:])
                    else:
                        matrix_t[i],matrix_p[i,j]=stats.ttest_rel(data[:,i],data1[:,j])
    im=ax.imshow(matrix_t,aspect='auto',extent=(0,len_data,0,len_data),cmap=cmap)
    significant_pixels=np.argwhere(matrix_p<p_level)
    
    for pix in significant_pixels:
        i=pix[0]
        j=pix[1]
        ax.axhline(len_data-i,j/len_data,(j+1)/len_data,color=color)
        ax.axhline(len_data-i-0.995,j/len_data,(j+1)/len_data,color=color)
        ax.axvline(j+0.01,(len_data-i-1)/len_data,(len_data-i)/len_data,color=color)
        ax.axvline(j+0.995,(len_data-i-1)/len_data,(len_data-i)/len_data,color=color)
        if annotate:
            ax.text(j+0.05,len_data-i-0.9,'%.1e'%matrix_p[i,j],fontsize=8)
    # print('Matrices')
    # print(matrix_t)
    # print(matrix_p)
    return im,ax

# data=np.array([[0,2,3,5,6,7,2,3,12,4,5,6,7,8],[16,25,34,85,22,28,32,45,15,26,28,29,30,20]])
# fig,ax=plt.subplots(1,1)
# im,ax=matrix_ttest(data,ax)
                
# datax=np.array([[2,3,4],[5,6,7],[8,9,10],[15,16,18]])
# data1x=np.array([[22,33,41],[15,16,17],[28,39,15],[17,26,38]])
# datay=np.array([[2,3,9],[12,6,7],[8,16,10],[34,16,18]])
# fig,ax=plt.subplots(1,2)
# ax[0]=univariable_std(datax,datay,ax[0],flagstdx=True,flagstdy=True,labels=['a','b','c'],relative_diference=False,meanzero_type='line')
# ax[1]=univariable_std(datax,datay,ax[1],flagstdx=True,flagstdy=True,cmap=plt.cm.Greens,labels=['a','b','c'],relative_diference=False)
# ax[1]=bivariable_std(datax,datay,ax[1],flagstdx=True,flagstdy=True,relative_diference=False)
# ax[1]=bivariable_std(data1x,datay,ax[1],flagstdx=True,flagstdy=True,relative_diference=False,cmap=plt.cm.Greens,marker='^')

# fig,ax=plt.subplots(1,2,subplot_kw=dict(projection='polar'))
# radii=np.random.rand(100)*np.pi
# ax[0]=circ_hist(radii,ax[0],bins=20) 
# ax[1]=circ_hist(radii,ax[1],bins=15)  
# ax[2]=circ_hist(radii,ax[2],bins=10)
# plt.show()
    