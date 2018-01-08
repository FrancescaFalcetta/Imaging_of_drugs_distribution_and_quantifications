#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
from scipy import misc
from scipy import ndimage
import scipy as sp
from PIL import Image
import matplotlib as mpl
import matplotlib.pyplot as plt
import os, struct
from numpy import *
import numpy as np
from skimage import filters
import argparse
import itertools
from pandas import *
import pandas as pd
import re
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import pandas.rpy.common as com
import csv
import time

ro.r('library(MALDIquant)') ##to call R function (i.e. removebaseline)


## Functions to read the Analyze 7.5 format
def readAnalyzeheader(filename):
   '''Function reading Analyze header file
   
   The header contains informations regardin gthe number of 
   points and other dataset characteristics.

   Args:
   filename : The name of the hdr file

   Value:
   nx : number of pixels in th x direction
   ny : number of pixel in the y direction
   xd : step in x 
   xy : step in y
   '''
   size = os.stat(filename).st_size
   if((size != 348) & (size != 384)):
       sys.exit("%s file is not Anlayze format.", filename)
   dimensions = np.zeros(8)
   pixdim = np.zeros(8)
   f = open(filename,'rb')
   f.seek(38)
   regular = f.read(1)
   if(regular != 'r'):
       sys.exit("wrong file format")
   f.seek(40)
   for i in range(8):
       dimensions[i] = struct.unpack('H',f.read(2))[0]
   nx = dimensions[2]  # x size of the image (number of columns)
   ny = dimensions[3]  # y size of the image (number of rows)
   f.seek(70)
   datatype = struct.unpack('H',f.read(2))[0]
   bitpix = struct.unpack('H',f.read(2))[0]
   if(datatype == 2 or datatype == 4 or datatype == 8 ):
       what = 'integer'
   elif (datatype == 16 or datatype == 32 or datatype == 64 ):
       what = 'double'
   else:
       what = 'raw'
   signed = (datatype == '2') 
   size = (bitpix/8)
   f.seek(76) 
   for i in range(8):
       pixdim[i] = struct.unpack('f',f.read(4))[0]
   xd = pixdim[1]
   yd = pixdim[2]
   return(int(nx),int(ny), xd, yd, signed, size, what) 

## read t2m
def readAnalyzet2m(filename):
    '''Function reading t2m mass file '''
    totalBytes = os.stat(filename).st_size
    Bytes = 4
    mass = np.zeros(totalBytes/Bytes)
    endian = 'f'
    with open(filename,'rb') as f:
        for i in xrange(len(mass)):
            mass[i] = struct.unpack(endian,f.read(Bytes))[0]
    return np.array(mass)  

 # exctact all spectra in the image   
def MSIcube(filename,t2m,hdr,samplename): 
    mass1, mass2 = min(t2m), max(t2m)
    np.savetxt(samplename + '_mass4cifre.sim',t2m,fmt='%-3.4f',delimiter=',') 
    np.savetxt(samplename + '_mass.sim',t2m,fmt='%-3.2f',delimiter=',') 
    idx = (t2m >= mass1) & (t2m <= mass2)
    Bytes, endian = 2, 'H' 
    Nmasses=sum(idx)
    cubematrix = np.zeros([Nmasses,hdr[0]* hdr[1]]) ## the matrix to store the image # 
    with open(filename,'rb') as f:
        for c in range(hdr[1] * hdr[0]): # all the pixels
            spcarr = []
            for i in range(t2m.size):
                if not idx[i]:
                    f.seek(Bytes,1)
                    continue 
                contenuto=f.read(Bytes)                    
                if contenuto=='':
                    spcarr=0
                    continue
                test = struct.unpack(endian,contenuto)[0]                
                spcarr.append(test)   
            cubematrix[:,c] = spcarr 
    return (cubematrix,Nmasses)

## read the image
def extraction_ions(idx,a,t2m,cubematrix,spcarr,area_sottr,mass1,mass2):   
    for i in range(t2m.size):
        if idx[i]:
            test=cubematrix[i,a]
            spcarr.append(test)
    if (area_sottr==[]): # calculation of the total area without background subtraction
        datapixel = sum(spcarr) 
    elif area_sottr=='range':  
        area_range=(spcarr[0]+spcarr[len(spcarr)-1])/2*len(spcarr)
        if area_range>sum(spcarr):
            datapixel=0
        else: ## subtraction from the total area, a trapezoidal background area   
            datapixel=sum(spcarr)-area_range 
    elif area_sottr=='min':
        area_min=np.min(spcarr)*len(spcarr) 
        if area_min>sum(spcarr):
            datapixel=0
        else:  ## subtract from the total area, a rectangular background area   
            datapixel=sum(spcarr)-area_min          
    return datapixel
    
def readAnalyzeImage(cubematrix,t2m,hdr,massrange = [],area_sottr = []):
    '''
    This function  interacts with readAnalyzet2m and readAnalyzeheader
    and extract 2d map showing the integrated intensity in the data
    '''
    datamatrix = np.zeros([hdr[0],hdr[1]]) ## the matrix to store the image
    a=0   
    for c,r in itertools.product(range(hdr[1]),range(hdr[0])):
            spcarr = []
            if (massrange == []):
                 mass1, mass2 = min(t2m), max(t2m)
                 idx=[]
                 datamatrix[r,c]=extraction_ions(idx,a,t2m,cubematrix,spcarr,area_sottr,mass1,mass2) 
            else:                     
                mass1, mass2 = massrange[a,0], massrange[a,1]
                if (max(massrange[a])!=0):     
                    idx = (t2m >= mass1) & (t2m <= mass2) 
                if (np.all(massrange[a])==0):
                    datamatrix[r,c]=0
                else:
                    datamatrix[r,c]=extraction_ions(idx,a,t2m,cubematrix,spcarr,area_sottr,mass1,mass2) 
            a=a+1                  
    return(datamatrix)

def getidmaxIntensity(filename,t2m):
    spectra = []
    f = open(filename,'rb')
    f.seek(t2m.size * int(20) * 2)
    for i in xrange(t2m.size):
        test = struct.unpack('H',f.read(2))[0]
        spectra.append(test)
    spectra = np.array(spectra)    
    nu =  np.where(spectra == max(spectra))[0][0]
    return(nu)

def grep(l,s):
    return [i for i in l if s in i]
    
def peakFound(cubematrix,t2m,hdr,loc_IS,delta_drug,delta_tiss,hws,hws_t,loc_tiss):

    devst=mean(t2m[1:t2m.size-1]-t2m[0:t2m.size-2]) 
    matrix_data=np.zeros([np.size(t2m),2]) 
    peakIS_matrix=np.zeros([hdr[0]*hdr[1],1])
    peakdrug_matrix=np.zeros([hdr[0]*hdr[1],1])
    peaktissue_matrix=np.zeros([hdr[0]*hdr[1],1])
    range_drugmatrix=np.zeros([hdr[0]*hdr[1],2])
    range_tissmatrix=np.zeros([hdr[0]*hdr[1],2])
    range_ISmatrix=np.zeros([hdr[0]*hdr[1],2])
    matrix_data[:,0]=t2m 
    
    for c in range(hdr[1] * hdr[0]): 
         matrix_data[:,1]=cubematrix[:,c]  # spettrum of single pixel "c"
         matrix_data_DF = pd.DataFrame(matrix_data,dtype='float')
         data = com.convert_to_r_matrix(matrix_data_DF)
         ro.globalenv['matrix_data_DF']=data 
         ro.globalenv['hws']=hws
         ro.globalenv['loc_IS']=loc_IS
         ro.globalenv['hws_t']=hws_t
         ro.globalenv['loc_tiss']=loc_tiss
         ro.globalenv['delta_drug']=delta_drug
         ro.globalenv['delta_tiss']=delta_tiss
         ro.globalenv['devst']=devst
         ro.r('spectra<-createMassSpectrum(mass=matrix_data_DF[,c(1)], intensity=matrix_data_DF[,c(2)])')
         ro.r('peaks<-detectPeaks(spectra, SNR=2,halfWindowSize=hws)') 
         ## identification of peak of internal standard
         ro.r('all_peaks_IS<-peaks@mass[(peaks@mass>=loc_IS[1])&(peaks@mass<=loc_IS[2])]')# all the peaks in this m/z range 
         ro.r('posISs<-match(all_peaks_IS,spectra@mass)')
         ro.r('NpeaksIS<-length(all_peaks_IS)')
         NpeaksIS=ro.r('NpeaksIS')[0]
         if NpeaksIS==0:
             continue
         if NpeaksIS>1:
             ro.r('info_peaksIS<-matrix(data = 0,nrow = NpeaksIS ,ncol =3 )')
             ro.r('info_peaksIS[,c(1)]<-posISs')
             ro.r('info_peaksIS[,c(2)]<-all_peaks_IS')
             ro.r('info_peaksIS[,c(3)]<-spectra@intensity[posISs]')
             ro.r('pos_peakIS_moda<-info_peaksIS[match(max(info_peaksIS[,c(3)]),info_peaksIS[,c(3)]),c(1)]')
         else:
             ro.r('pos_peakIS_moda<-posISs')
         ro.r('pos_peakIS_moda')[0]    
         ro.r('matrix_lav<-matrix(data = 0,nrow = hws*2+1 ,ncol =4 )')
         ro.r('matrix_lav[,c(1)]<-spectra@mass[(pos_peakIS_moda-hws):(pos_peakIS_moda+hws)]')
         ro.r('matrix_lav[,c(2)]<-spectra@intensity[(pos_peakIS_moda-hws):(pos_peakIS_moda+hws)]')
         ro.r('matrix_lav[,c(3)]<-(pos_peakIS_moda-hws):(pos_peakIS_moda+hws)')
         ro.r('peakIS_mean<-sum(matrix_lav[,c(1)]*matrix_lav[,c(2)])/sum(matrix_lav[,c(2)])')
         ro.r('matrix_lav[,c(4)]<-abs(matrix_lav[,c(1)]-peakIS_mean)')
         pos_peakIS=ro.r('pos_peakIS<-matrix_lav[match(min(matrix_lav[,c(4)]),matrix_lav[,c(4)]),c(3)]')[0]
         peaks_IS=ro.r('peakIS_mean')
         if np.isnan(peaks_IS):
             continue
         if NpeaksIS>0:
             ro.r('range_IS<-c(spectra@mass[pos_peakIS-hws], spectra@mass[pos_peakIS+hws])')
             rim=ro.r('range_IS')
             range_ISmatrix[c,0]=rim[0]
             range_ISmatrix[c,1]=rim[1]
             peakIS_matrix[c]=peaks_IS[0]
             ro.r('all_peaks_drug<-spectra@mass[(spectra@mass>=peakIS_mean-delta_drug-devst)&(spectra@mass<=peakIS_mean-delta_drug+devst)]')
             ro.r('Npeaks_drug<-length(all_peaks_drug)')
             all_peaks_drug=ro.r('all_peaks_drug')
             Npeaks_drug=ro.r('Npeaks_drug')[0]
             if Npeaks_drug>1:
                 peaks_drug=0
                 p=0
                 min_diff_peaks_drug=ro.r('min(abs(peakIS_mean-delta_drug-all_peaks_drug))')[0]
                 while peaks_drug<=0:
                     if abs(peaks_IS[0]-delta_drug-all_peaks_drug[p])==min_diff_peaks_drug:
                            peaks_drug=all_peaks_drug[p]  
                     p=p+1
             else:
                 peaks_drug=all_peaks_drug
             peakdrug_matrix[c]= peaks_drug   
             ro.globalenv['peaks_drug']=peaks_drug
             ro.r('pos_peakdrug<-match(peaks_drug,spectra@mass)')
             ro.r('range_drug<-c(spectra@mass[pos_peakdrug-hws], spectra@mass[pos_peakdrug+hws])')
             rdm=ro.r('range_drug')
             range_drugmatrix[c,0]=rdm[0]
             range_drugmatrix[c,1]=rdm[1]
             # identification of peak of tissue
             ro.r('all_peaks_tiss<-peaks@mass[(peaks@mass>=loc_tiss[1])&(peaks@mass<=loc_tiss[2])]')# tutti i picchi nell'intervallo cercato
             ro.r('pos_tiss<-match(all_peaks_tiss,spectra@mass)')
             ro.r('Npeaks_tiss<-length(all_peaks_tiss)')
             NpeaksTiss=ro.r('Npeaks_tiss')[0]
             if NpeaksTiss==0:
                 continue
             if NpeaksTiss>1:
                 ro.r('info_peaksTiss<-matrix(data = 0,nrow = Npeaks_tiss ,ncol =3 )')
                 ro.r('info_peaksTiss[,c(1)]<-pos_tiss')
                 ro.r('info_peaksTiss[,c(2)]<-all_peaks_tiss')
                 ro.r('info_peaksTiss[,c(3)]<-spectra@intensity[pos_tiss]')
                 ro.r('peaks_tiss<-info_peaksTiss[match(max(info_peaksTiss[,c(3)]),info_peaksTiss[,c(3)]),c(2)]')
                 ro.r('pos_peakTiss_moda<-info_peaksTiss[match(max(info_peaksTiss[,c(3)]),info_peaksTiss[,c(3)]),c(1)]')
             else:
                 ro.r('pos_peakTiss_moda<-pos_tiss')
             ro.r('pos_peakTiss_moda')[0]    
             ro.r('matrix_lav<-matrix(data = 0,nrow = hws_t*2+1 ,ncol =4 )')
             ro.r('matrix_lav[,c(1)]<-spectra@mass[(pos_peakTiss_moda-hws_t):(pos_peakTiss_moda+hws_t)]')
             ro.r('matrix_lav[,c(2)]<-spectra@intensity[(pos_peakTiss_moda-hws_t):(pos_peakTiss_moda+hws_t)]')
             ro.r('matrix_lav[,c(3)]<-(pos_peakTiss_moda-hws_t):(pos_peakTiss_moda+hws_t)')
             ro.r('peaktiss_mean<-sum(matrix_lav[,c(1)]*matrix_lav[,c(2)])/sum(matrix_lav[,c(2)])')
             ro.r('matrix_lav[,c(4)]<-abs(matrix_lav[,c(1)]-peaktiss_mean)')    
             pos_peakTiss=ro.r('pos_peaktiss<-matrix_lav[match(min(matrix_lav[,c(4)]),matrix_lav[,c(4)]),c(3)]')[0]
             peaks_tiss=ro.r('peaktiss_mean')
             if np.isnan(peaks_tiss):
                 continue
             if NpeaksTiss>0:
                 ro.r('range_tiss<-c(spectra@mass[pos_peaktiss-hws_t], spectra@mass[pos_peaktiss+hws_t])')
                 rtm=ro.r('range_tiss')
                 range_tissmatrix[c,0]=rtm[0]
                 range_tissmatrix[c,1]=rtm[1]  
                 peaktissue_matrix[c]=peaks_tiss[0]
    return (peakIS_matrix,peakdrug_matrix,peaktissue_matrix,range_ISmatrix, range_drugmatrix, range_tissmatrix)

## batch process stuff ...
def processMSIBatch(mfiltrad = 3,val_in=[],loc_IS=[],delta_drug=[],delta_tiss=[],hws=5,scale_factors=1.0,area_sottr = [],r_p=4,c_p=4, hws_t=[],loc_tiss=[]):
    '''
    Process in batch mode all the analyze files contained in a FOlder
    '''      
    tic= time.clock()
    dircontent = os.listdir('.')  ## directory 
    hdrs = grep(dircontent,'.hdr')
    val_ctrl=[]
    for f in hdrs:
        MSImatrixCube=[]
        MSImatrix_drug=[]
        MSImatrix_tissue=[]
        MSImatrix_std=[]
        mask_drug=[]
        mask_tissue=[]
        flattened_list=[]     
        myheader = readAnalyzeheader(f)
        sname = f[:-4]      ## get rid of the extension
        print('PROCESSING %s : ' % (sname))      
        mass = readAnalyzet2m(sname + '.t2m')
        ## extract and save the matrix containing all spectra
        MSImatrixCube,Nmasses=MSIcube(sname + '.img',mass,myheader,sname)
        np.savetxt(sname +'_'+ 'MSImatrixCube.sim',MSImatrixCube,fmt='%-3.2f',delimiter=',')   
        if hws_t==[]:
            hws_t=hws
            
######## peaks identification    
        peakstd_matrix,peakdrug_matrix,peaktissue_matrix,massrange_std_matrix,massrange_drug_matrix,massrange_tiss_matrix=peakFound(MSImatrixCube,mass,myheader,loc_IS,delta_drug,delta_tiss,hws,hws_t,loc_tiss)
        massrange=[]
        massrange_tissue=[]
        massrange_std=[]
        massrange=massrange_drug_matrix
        massrange_tissue=massrange_tiss_matrix
        massrange_std=massrange_std_matrix
        np.savetxt(sname +'_'+ 'peakstd.sim',peakstd_matrix,fmt='%-3.2f',delimiter=',')   
        np.savetxt(sname +'_'+ 'peakdrug.sim',peakdrug_matrix,fmt='%-3.2f',delimiter=',')     
        np.savetxt(sname +'_'+ 'peaktissue.sim',peaktissue_matrix,fmt='%-3.2f',delimiter=',')      
        np.savetxt(sname +'_'+ 'massrange_IS.sim',massrange_std_matrix,fmt='%-3.2f',delimiter=',')  
        np.savetxt(sname +'_'+ 'massrange_drug.sim',massrange_drug_matrix,fmt='%-3.2f',delimiter=',')  
        np.savetxt(sname +'_'+ 'massrange_tiss.sim',massrange_tiss_matrix,fmt='%-3.2f',delimiter=',')      
        
#prodecure with ion of tissue
        if np.any(massrange_tissue): 
######## peaks integration
               MSImatrix_drug = readAnalyzeImage(MSImatrixCube,mass,myheader,massrange,area_sottr) 
               np.savetxt(sname +'_'+ 'MSImatrix_drug.sim',MSImatrix_drug,fmt='%-3.2f',delimiter=',')   
               MSImatrix_tissue = readAnalyzeImage(MSImatrixCube,mass,myheader,massrange_tissue,area_sottr)               
               np.savetxt(sname +'_'+ 'MSImatrix_tissue.sim',MSImatrix_tissue,fmt='%-3.2f',delimiter=',')   

############## Normalization Dividing drug and tissue image with standard              
               if np.any(massrange_std):
                   MSImatrix_std = readAnalyzeImage(MSImatrixCube,mass,myheader,massrange_std,area_sottr) 
                   np.savetxt(sname +'_'+ 'MSImatrix_std.sim',MSImatrix_std,fmt='%-3.2f',delimiter=',')   
                   MSImatrix_std[MSImatrix_std==0]=1
                   MSImatrix_drug = MSImatrix_drug / MSImatrix_std
                   MSImatrix_tissue = MSImatrix_tissue / MSImatrix_std
                             
############## Tissue Isolation  
               print('<--- Extracting and saving tissue object -->')     
               mask_tissue_bin=[]
               MSImatrix_tissue = ndimage.median_filter(MSImatrix_tissue,mfiltrad) #noise reduction
               nr = MSImatrix_tissue.shape[0]-1
               nc = MSImatrix_tissue.shape[1]-1
               ul=(MSImatrix_tissue[0:r_p,0:c_p])
               ur=(MSImatrix_tissue[0:r_p,(nc-c_p+1):nc+1])
               dl=(MSImatrix_tissue[(nr-r_p+1):nr+1,0:c_p])
               dr=(MSImatrix_tissue[(nr-r_p+1):nr+1,(nc-c_p+1):nc+1])
               nu=[]    
               nu.append(ul)
               nu.append(ur)  
               nu.append(dr)
               nu.append(dl)   
               val_piastra=np.percentile(nu, 99) 
               mask_tissue = MSImatrix_tissue > val_piastra
               mask_tissue_bin=mask_tissue.astype(int) #modificato per produrre maschera binaria
               mask_tissue = mask_tissue.astype(int)
               np.savetxt(sname + 'mask_tissue_bin.msk',mask_tissue_bin,fmt ='%-3.2f',delimiter=',')  
               
############## Create drug mask   
               print('<--- Exracting and Saving Drug Image -->')
               mask_drug = MSImatrix_drug * mask_tissue
               mask_tissue = ndimage.median_filter(mask_drug,mfiltrad)
               mask_drug = ndimage.median_filter(mask_drug,mfiltrad) 
               if val_in ==[]:
                   patterns=['CTRL','ctrl','Ctrl']
                   a=np.zeros(3,dtype=bool)
                   i=0
                   for pattern in patterns:                      
                       a[i]=re.search(pattern, sname)
                       i=i+1
                   if np.any(a):
                         val_ctrl = np.percentile(mask_drug, 95) #                  
                         print('val ctrl 95 percentile %-3.2f : ' % (val_ctrl))
                   if val_ctrl==[]:
                       val = filters.threshold_otsu(mask_drug)
                       print('val_otsu %-3.2f : ' % (val))
                   else:
                       val=val_ctrl
               else:
                   val=val_in
               mask_tissue = mask_drug> val
          
               np.savetxt(sname + '.sim',mask_drug,fmt='%-3.2f',delimiter=',')
               np.savetxt(sname + '_drug.msk',mask_tissue,fmt ='%-3.2f',delimiter=',')
               fig = plt.figure()
               fig.subplots_adjust(left=None, bottom=None, right=None, top=None,wspace=None, hspace=0.3)
               ax1 = plt.subplot(311)
               plt.title(sname)
               plt.imshow(mask_drug,interpolation='None', clim=(0.0, scale_factors))
               plt.colorbar()
               ax2 = plt.subplot(312)
               plt.title('drug mask')
               plt.imshow(mask_tissue,interpolation='None',cmap='Greys')
               plt.colorbar()
               ax3 = plt.subplot(313)
               plt.title('tissue mask')
               plt.imshow(mask_tissue_bin,interpolation='None',cmap='Greys')
               plt.colorbar()
               plt.show()
               plt.savefig(sname + '.jpg',bbox_inches='tight')
               plt.close()
               toc=time.clock()
               elapsed_time=toc-tic       
               print('elapsed time %-3.2f : '%(elapsed_time))
        else: #without ion of tissue  
               print('<--- Exracting and Saving Drug Image -->')
               
               ######## peaks integration
               MSImatrix = readAnalyzeImage(sname + '.img',mass,myheader,massrange,area_sottr) 
                         
               if np.any(massrange_std):
                   MSImatrix_std = readAnalyzeImage(sname + '.img',mass,myheader,massrange_std,area_sottr) 
                   ##### normalization
                   MSImatrix_std[MSImatrix_std==0]=1
                   MSImatrix = MSImatrix / MSImatrix_std
                   ### noise reduction
                   MSImatrix = ndimage.median_filter(MSImatrix,mfiltrad)

############## Create drug mask
                   print('<--- Exracting and Saving Drug Mask -->')
                   if val_in ==[]:
                       patterns=['CTRL','ctrl','Ctrl']
                       val_otsu = filters.threshold_otsu(mask_drug)
                       a=np.zeros(3,dtype=bool)
                       i=0
                       for pattern in patterns:                           
                           a[i]=re.search(pattern, sname)
                           i=i+1
                       if np.any(a):
                            val_ctrl = np.percentile(mask_drug, 95) #                  
                            print('val ctrl 95 percentile %-3.2f: ' % (val_ctrl))
                       if val_ctrl==[]:
                           val = val_otsu
                           print('val_otsu %-3.2f : ' % (val))
                       else:
                           val=val_ctrl
                           print('val ctrl %-3.2f vs val_Otsu %-3.2f : ' % (val_ctrl, val_otsu))
                   else:
                       val=val_in
                   mask_tissue = mask_drug> val
                   mask_drug = np.sqrt(MSImatrix) > val
                   mask_drug = mask_drug.astype(int)        
                   np.savetxt(sname + '.sim',MSImatrix,fmt='%-3.2f',delimiter=',')
                   np.savetxt(sname + '_drug.msk',mask_drug,fmt ='%-3.2f',delimiter=',')
                   fig = plt.figure(1)
                   ax1 = plt.subplot(211)
                   plt.title(sname)
                   plt.imshow(MSImatrix,interpolation='None')
                   plt.colorbar()
                   ax2 = plt.subplot(212)
                   plt.title('drug mask')
                   plt.imshow(mask_drug,interpolation='None',cmap='Greys')
                   plt.colorbar()
                   plt.savefig(sname + '.jpg',bbox_inches='tight')
                   plt.close()

def main():    
    parser = argparse.ArgumentParser(description="Create an extracted ion image and its corresponding binary mask")
    parser.add_argument('-fr',dest = "mfiltrad",type = int, default=3,help = "Radius of the median filter")
    parser.add_argument('-th_mask',dest = "val_in",type = float,default=[],help = "threshold value for the drug mask, given as input or automatically calculated from control sample")  
    parser.add_argument('-locIS', dest = "loc_IS",type = float, nargs ='+', default=[],help = "m/z range where to find the internal standard (IS) peak")  
    parser.add_argument('-d_drug',dest = "delta_drug",type = float,default=[],help = "m/z difference between IS and drug peaks") 
    parser.add_argument('-locTiss',dest = "loc_tiss",type = float, nargs ='+',default=[],help =  "m/z range where to find the tissue peak")     
    parser.add_argument('-hws',dest = "hws",type = float,default=5 ,help = "half -width size of IS, drug and tissue peaks . Default value 5") 
    parser.add_argument('-hws_t',dest = "hws_t",type = float,default=[] ,help = "half -width size of tissue, peak used only when diefferent from -hws ") 
    parser.add_argument('-fs',dest = "scale_factors",type = float,default=1.0,help = "Maximum intensity in the color map. (Default fs=1.0)")    
    parser.add_argument('-a',dest = "area_sottr",type = str,default=[],help = "modality to calculate the peak area. Options: a=[]:No background subtraction.-a range:trapezoidal background subtraction.-a min:rectangular background subtraction.(Default a=[])")    
    parser.add_argument('-rp',dest="r_p",type=int, default=4, help="size (number of rows) of slice's edge used to estimate the plate intensity of the tissue ion peak")
    parser.add_argument('-cp',dest="c_p",type=int, default=4, help="size (number of columns) of slice's edge used to estimate the plate intensity of the tissue ion peak")
    args = parser.parse_args()
    processMSIBatch(mfiltrad = args.mfiltrad, val_in=args.val_in,loc_IS=args.loc_IS,delta_drug=args.delta_drug,loc_tiss=args.loc_tiss,hws=args.hws,hws_t=args.hws_t,scale_factors=args.scale_factors,area_sottr=args.area_sottr,r_p=args.r_p,c_p=args.c_p)
                    
if __name__ == '__main__':
   main()