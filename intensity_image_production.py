# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 09:13:26 2016

@author: FFalcetta
"""

from scipy import ndimage
import matplotlib.pyplot as plt
import os
from numpy import *
import numpy as np
from skimage import filters
import argparse
import re
import csv
import statsmodels.api as sm
from PIL import Image
from io import BytesIO

def grep(l,s):
    return [i for i in l if s in i]
    
def factoryscale(selionimg_ord):        
    Int_vals=[]
    for sf in selionimg_ord[1:]:
        sname = sf[:-4]
        MSImatrix_drug = np.genfromtxt(sname+'_MSImatrix_drug.sim',dtype=float,delimiter=',')
        MSImatrix_std=np.genfromtxt(sname+'_MSImatrix_std.sim',dtype=float,delimiter=',')
        MSImatrix_std[MSImatrix_std==0]=1
        MSImatrix_drug = MSImatrix_drug / MSImatrix_std
        int_val=np.mean(MSImatrix_drug)

        Int_vals.append(int_val)
        scale_factors=np.mean(Int_vals)
    return scale_factors

def smooth_medianFilter(sname,MSImatrix_drug,MSImatrix_std,MSImatrix_tissue, mfiltrad,val_in,scale_factors,val_ctrl,r_p,c_p):
### production of binary mask.    
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
    val_piastra=np.percentile(nu,99)
    mask_tissue_bin= MSImatrix_tissue > val_piastra
    mask_tissue_bin=mask_tissue_bin.astype(int)
    np.savetxt(sname + 'mask_tissue_bin.msk',mask_tissue_bin,fmt ='%-3.2f',delimiter=',')       
    mask_drug = MSImatrix_drug * mask_tissue_bin
    mask_drug = ndimage.median_filter(mask_drug,mfiltrad) 
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
        if val_ctrl==[]:
            val = val_otsu
            print('Quantification matrix: val_otsu %-3.2f ' % (val))
        else:
            val=val_ctrl
            print('Quantification matrix: val ctrl %-3.2f vs val_Otsu %-3.2f ' % (val_ctrl, val_otsu))
    else:
        val=val_in
    mask_tissue = mask_drug> val
    mask_tissue = mask_tissue.astype(int)
    mask_drug=mask_drug*mask_tissue  
    np.savetxt(sname + '.sim',mask_drug,fmt='%-3.2f',delimiter=',')
    np.savetxt(sname + '_drug.msk',mask_tissue,fmt ='%-3.2f',delimiter=',')
    return(mask_drug,mask_tissue,mask_tissue_bin,val_ctrl)

def figures_production(sname,mask_drug,mask_tissue,mask_tissue_bin,scale_factors,color_map):   
    pixel_tot_tissue=float64(sum(mask_tissue_bin))
    black_pixel=float64(sum(mask_tissue*mask_tissue_bin))
    percentage=black_pixel/pixel_tot_tissue*100
    intensity_tot=float64(sum(mask_drug*mask_tissue_bin))
    drug_IS=intensity_tot/pixel_tot_tissue
    np.savetxt(sname + '_drug.sim',mask_drug,fmt ='%-3.2f',delimiter=',')
    fig=plt.figure(dpi=300)
    if color_map:
        plt.imshow(mask_drug,aspect='equal',interpolation='None',cmap=color_map,clim=(0.0, scale_factors),origin='upper')
    else:
        plt.imshow(mask_drug,aspect='equal',interpolation='None', clim=(0.0, scale_factors),origin='upper')
    plt.colorbar()
    plt.axis('off')
    png1 = BytesIO()
    fig.savefig(png1, format='png',bbox_inches='tight', transparent='true')
    png2 = Image.open(png1)    
    png2.save(sname + 'drug_intensity.tiff')
    plt.close()

    fig=plt.figure(dpi=300)
    plt.imshow(mask_tissue,interpolation='None',cmap='Greys')
    plt.axis('off')
    png1 = BytesIO()
    fig.savefig(png1, format='png',bbox_inches='tight', transparent='true')
    png2 = Image.open(png1)    
    png2.save(sname + 'drug_mask.tiff')
    plt.close()
    
    fig=plt.figure(dpi=300)
    plt.imshow(mask_tissue_bin,interpolation='None',cmap='Greys')
    plt.axis('off')
    plt.show()
    png1 = BytesIO()
    fig.savefig(png1, format='png',bbox_inches='tight', transparent='true')
    png2 = Image.open(png1)    
    png2.save(sname + 'tissue_mask.tiff')
    plt.close()
    return(pixel_tot_tissue,black_pixel,percentage,drug_IS)
    
def process(val_in=[],scale_factors=[],r_p=4,c_p=4,color_map=[]):
    val_ctrl=[]
    all_pixel_tot_tissue=['pixel_tot_tissue']
    all_black_pixel=['black_pixel']
    all_percentage=['%']
    all_intensity_tot=['total intensity']
    vett_name=[' ']
    mfiltrad=3
    dircontent = os.listdir('.')

    selionimg = grep(dircontent,'.jpg') 
    
    patterns_ctrl=['CTRL','ctrl','Ctrl']
    n_ctrl=0
   
    for f in selionimg:
        a_c=np.zeros(3,dtype=bool)
        i_c=0
        sname = f[:-4] 
        for pattern_ctrl in patterns_ctrl:           
            a_c[i_c]=re.search(pattern_ctrl, sname)
            i_c=i_c+1
        if np.any(a_c):
            pos_ctrl=n_ctrl
        n_ctrl=n_ctrl+1    
    if pos_ctrl==0: 
        selionimg_ord=selionimg
    elif pos_ctrl==len(selionimg)-1:
        selionimg_ord=[selionimg[pos_ctrl]]+selionimg[:pos_ctrl]
    else:
        selionimg_ord=[selionimg[pos_ctrl]]+selionimg[:pos_ctrl]+selionimg[pos_ctrl+1:]      
    print selionimg_ord 
    
    if scale_factors==[]:
        scale_factors=factoryscale(selionimg_ord)
        
    for f_ord in selionimg_ord:
        sname = f_ord[:-4]
        a_c=np.zeros(3,dtype=bool)
        i_c=0
        print('PROCESSING %s : ' % (sname))

        MSImatrix_drug = np.genfromtxt(sname+'_MSImatrix_drug.sim',dtype=float,delimiter=',')
        MSImatrix_tissue=np.genfromtxt(sname+'_MSImatrix_tissue.sim',dtype=float,delimiter=',')
        MSImatrix_std=np.genfromtxt(sname+'_MSImatrix_std.sim',dtype=float,delimiter=',')
        mask_drug,mask_tissue,mask_tissue_bin,val_ctrl=smooth_medianFilter(sname,MSImatrix_drug,MSImatrix_std,MSImatrix_tissue, mfiltrad,val_in,scale_factors,val_ctrl,r_p,c_p)

        for pattern_ctrl in patterns_ctrl:           
            a_c[i_c]=re.search(pattern_ctrl, sname)
            i_c=i_c+1
        if np.any(a_c):
            sname_CTRL=sname
            mask_drug_CTRL=mask_drug
            mask_tissue_CTRL=mask_tissue
            mask_tissue_bin_CTRL=mask_tissue_bin            
            vett_name=vett_name+[sname_CTRL]
            pixel_tot_tissue,black_pixel,percentage,intensity_tot=figures_production(sname_CTRL,mask_drug_CTRL,mask_tissue_CTRL,mask_tissue_bin_CTRL,scale_factors,color_map)
        else:
            vett_name=vett_name+[sname]
            pixel_tot_tissue,black_pixel,percentage,intensity_tot=figures_production(sname,mask_drug,mask_tissue,mask_tissue_bin,scale_factors,color_map)
        all_pixel_tot_tissue=all_pixel_tot_tissue+[pixel_tot_tissue]
        all_black_pixel=all_black_pixel+[black_pixel]
        all_percentage=all_percentage+[percentage]
        all_intensity_tot=all_intensity_tot+[intensity_tot]
    results=[(vett_name),(all_pixel_tot_tissue),(all_black_pixel), (all_percentage),(all_intensity_tot)]
    with open('results.csv','wb')as csvfile:
        writer=csv.writer(csvfile)
        writer.writerows(results)

def main():    
    parser = argparse.ArgumentParser(description="drug intensity image production.")
    parser.add_argument('-th_mask',dest = "val_in",type = float,default=[],help = "threshold to obtain the drug_mask. Default [] Otsu’s method")  
    parser.add_argument('-rp',dest="r_p",type=float, default=4, help="size (number of rows) of slice's edge used to estimate the plate intensity of the tissue ion peak")
    parser.add_argument('-cp',dest="c_p",type=float, default=4, help="size (number of columns) of slice's edge used to estimate the plate intensity of the tissue ion peak")
    parser.add_argument('-co',dest="color_map",type = str,default=[],help = "choice of colormap type.") 
    parser.add_argument('-fs',dest = "scale_factors",type = float,default=[],help = "Maximum intensity in the color map. (Default fs=1.0)")
    args = parser.parse_args()
    process(val_in=args.val_in,scale_factors=args.scale_factors,r_p=args.r_p,c_p=args.c_p,color_map=args.color_map)
    
if __name__ == '__main__':
   main()    