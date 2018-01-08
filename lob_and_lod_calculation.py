# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 09:13:26 2016

@author: FFalcetta
"""

from scipy import ndimage
import os
from numpy import *
import numpy as np
from skimage import filters
from skimage import img_as_ubyte
import argparse
import re
import statsmodels.api as sm
from PIL import Image
from io import BytesIO

def grep(l,s):
    return [i for i in l if s in i]   
    
def lob_lod(sname,mask_drug,mask_tissue_bin,coord_spot,conc_drug):
    nspot=len(coord_spot)/4
    if nspot==len(conc_drug):
        Int_spot=np.zeros(nspot)
        dim_spot=np.zeros(nspot)
        conc_spot_im=np.zeros(nspot) 
        conc_spot=np.zeros(nspot) 
        perc5_spot=np.zeros(nspot)
        mean_Int_spot=np.zeros(nspot)
        ###### identification of balnk spot and LOB calculation
        for z in range(nspot):
            if conc_drug[z]==0:
                 ri=coord_spot[z*4]
                 rf=coord_spot[z*4+1]
                 ci=coord_spot[z*4+2]
                 cf=coord_spot[z*4+3]
                 lob=np.percentile(mask_drug[ri:rf+1,ci:cf+1],95)
                 pos0=z
        if pos0==0:
            nspot_f=nspot
            p=1
        if pos0==nspot-1:
            nspot_f=nspot-1
            p=0
        ##### identification other spots and LOD calculation    
        for s in range(p,nspot_f):
            ri=coord_spot[(s)*4]
            rf=coord_spot[(s)*4+1]
            ci=coord_spot[(s)*4+2]
            cf=coord_spot[(s)*4+3]
            Int_spot[s]=np.sum(mask_drug[ri:rf+1,ci:cf+1])
            dim_spot[s]=np.sum(mask_tissue_bin[ri:rf+1,ci:cf+1])
            conc_spot_im[s]=Int_spot[s]/dim_spot[s] #y
            conc_spot[s]=conc_drug[s]/dim_spot[s] #x
            perc5_spot[s]=np.percentile(mask_drug[ri:rf+1,ci:cf+1],5)
            mean_Int_spot[s]=np.mean(mask_drug[ri:rf+1,ci:cf+1])  

        if pos0==0:
            diff_lob=perc5_spot[1:]-lob
        if pos0==nspot-1:
            diff_lob=perc5_spot[:nspot_f]-lob
        dif_lob_abs=np.abs(diff_lob)
        min_diff=np.min(dif_lob_abs)
        for h in range(len(diff_lob)):
            if diff_lob[h]==min_diff:
                pos_lod=h
        if pos0==0:
            pos_lod=pos_lod+1
        if perc5_spot[pos_lod]==lob:
            lod=mean_Int_spot[pos_lod]
        else:
            lod=mean_Int_spot[pos_lod]-diff_lob[pos_lod]
        conc_spot=conc_drug/(dim_spot*0.01) #expected concentration x /mm2

    return lob,lod

def smooth_medianFilter(sname,MSImatrix_drug,MSImatrix_std,MSImatrix_tissue, mfiltrad,r_p,c_p):
    mask_tissue_bin=[]
    MSImatrix_std[MSImatrix_std==0]=1
    MSImatrix_drug = MSImatrix_drug / MSImatrix_std
    MSImatrix_tissue = np.sqrt(MSImatrix_tissue) / (np.sqrt(MSImatrix_std)+1)
    MSImatrix_tissue = ndimage.median_filter(MSImatrix_tissue,mfiltrad)
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
    val_piastra=np.percentile(nu, 95) 
    mask_tissue_bin= MSImatrix_tissue > val_piastra
    mask_tissue_bin=mask_tissue_bin.astype(int)
    np.savetxt(sname + 'mask_tissue_bin.msk',mask_tissue_bin,fmt ='%-3.2f',delimiter=',')     
    mask_drug = MSImatrix_drug * mask_tissue_bin
    mask_drug = ndimage.median_filter(mask_drug,mfiltrad) 
    np.savetxt(sname + '.sim',mask_drug,fmt='%-3.2f',delimiter=',')

    return(mask_drug,mask_tissue_bin)


def process(coord_spot=[],conc_drug=[],r_p=4,c_p=4):
    mfiltrad=3
    dircontent = os.listdir('.')
    selionimg = grep(dircontent,'.jpg')
    patterns_retta=['CALIBRATION','calibration','Calibration']
    for f in selionimg:
        a_r=np.zeros(9,dtype=bool)
        i_r=0
        sname = f[:-4] 
        for pattern_retta in patterns_retta:                      
            a_r[i_r]=re.search(pattern_retta, sname)
            i_r=i_r+1
        if np.any(a_r):            
            print('PROCESSING %s : ' % (sname))
            MSImatrix_drug = np.genfromtxt(sname+'_MSImatrix_drug.sim',dtype=float,delimiter=',')
            MSImatrix_tissue=np.genfromtxt(sname+'_MSImatrix_tissue.sim',dtype=float,delimiter=',')
            MSImatrix_std=np.genfromtxt(sname+'_MSImatrix_std.sim',dtype=float,delimiter=',')
            ###preprocessing steps 
            mask_drug,mask_tissue_bin=smooth_medianFilter(sname,MSImatrix_drug,MSImatrix_std,MSImatrix_tissue, mfiltrad,r_p,c_p)
            ###lob and lod calculation
            lob,lod=lob_lod(sname,mask_drug,mask_tissue_bin,coord_spot,conc_drug)
            print('Lob:  %-3.4f ' % (lob))
            print('Lod:  %-3.4f ' % (lod))

def main():    
    parser = argparse.ArgumentParser(description="script to calculate Lob and Lod value into calibration curve sample")
    parser.add_argument('-sp',dest="coord_spot",type=int, nargs ='+', default=[], help="localization of spots of the calibration curve into the tissue, i.e. the row and column numbers of the corners (four numbers ) for each spot")
    parser.add_argument('-c',dest="conc_drug",type=float, nargs ='+', default=[], help="drug concentration of calibration spots")
    parser.add_argument('-rp',dest="r_p",type=float, default=4, help="size (number of rows) of slice's edge used to estimate the plate intensity of the tissue ion peak")
    parser.add_argument('-cp',dest="c_p",type=float, default=4, help="size (number of columns) of slice's edge used to estimate the plate intensity of the tissue ion peak")
    args = parser.parse_args()
    process(coord_spot=args.coord_spot,conc_drug=args.conc_drug,r_p=args.r_p,c_p=args.c_p)
    
if __name__ == '__main__':
   main()    