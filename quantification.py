# -*- coding: utf-8 -*-
"""
Created on Wed Nov 30 09:13:26 2016

@author: FFalcetta
"""

from scipy import ndimage
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.cbook as cbook
#from matplotlib_scalebar.scalebar import ScaleBar
import matplotlib.colors as colors
import os, struct
from numpy import *
import numpy as np
from skimage import filters
from skimage import img_as_ubyte
import argparse
import re
import csv
import statsmodels.api as sm
import PIL
from PIL import Image
from io import BytesIO

def grep(l,s):
    return [i for i in l if s in i]   
  
def figures_production(sname,mask_drug_conc,mask_tissue,mask_tissue_bin,scale_factors,color_map,voxel,lod_conc,unit_F):   

    fig=plt.figure(dpi=300)
    if color_map:
        plt.imshow(mask_drug_conc,aspect='equal',interpolation='None', cmap=color_map, vmin=lod_conc, vmax=scale_factors,origin='upper')
    else:
        plt.imshow(mask_drug_conc,aspect='equal',interpolation='None', vmin=lod_conc, vmax=scale_factors,origin='upper')
    cbar=plt.colorbar()
    cbar.ax.set_xlabel(unit_F)
    plt.show()    
    plt.axis('off')
    png1 = BytesIO()
    fig.savefig(png1, format='png',bbox_inches='tight', transparent='true')
    png2 = Image.open(png1)
    png2.save(sname + 'drug_concentration.tiff')
    png1.close()
    plt.close()     
    fig=plt.figure()
    plt.title(sname+'drug mask')
    plt.imshow(mask_tissue,interpolation='None',cmap='Greys')
    plt.axis('off')
    png1 = BytesIO()
    fig.savefig(png1, format='png')
    png2 = Image.open(png1)
    png2.save(sname + 'drug_mask.tiff')
    plt.close()
    png1.close()
    fig=plt.figure()
    plt.title(sname+'tissue mask')
    plt.imshow(mask_tissue_bin,interpolation='None',cmap='Greys')
    plt.axis('off')
    png1 = BytesIO()
    fig.savefig(png1, format='png')
    png2 = Image.open(png1)
    png2.save(sname + 'tissue_mask.tiff')
    plt.close()
    png1.close()   
    return()
   
def calibration_curve(sname,mask_drug,mask_tissue,coord_spot,conc_drug,weight_type,mw_drug,scale_factors,color_map,unit_F):
    nspot=len(coord_spot)/4
    if nspot==len(conc_drug):
        Int_spot=np.zeros(nspot)
        dim_spot=np.zeros(nspot)
        conc_spot_im=np.zeros(nspot) #y
        conc_spot=np.zeros(nspot) #x
        for s in range(nspot):
            ri=coord_spot[s*4]
            rf=coord_spot[s*4+1]
            ci=coord_spot[s*4+2]
            cf=coord_spot[s*4+3]
            Int_spot[s]=np.sum(mask_drug[ri:rf+1,ci:cf+1])
            dim_spot[s]=np.sum(mask_tissue[ri:rf+1,ci:cf+1])
            conc_spot_im[s]=Int_spot[s]/dim_spot[s] #y
            conc_spot[s]=conc_drug[s]/(dim_spot[s]*0.01) #x conc attesa/mm2
        if weight_type:
            if weight_type==1:
                weight=1/conc_spot_im
                weight[weight ==inf]=0
                weight_str="1/y^2"
            elif weight_type==2:
                weight=1/conc_spot
                weight[weight ==inf]=0
                weight_str="1/x^2"
            coef=np.polyfit(conc_spot, conc_spot_im,1,w=weight)
        else:
            coef=np.polyfit(conc_spot, conc_spot_im,1)
            weight_str="none"
        m=coef[0]
        intercept=coef[1]

        yth=m*conc_spot + intercept # mm2
        fig=plt.figure(dpi=300)
        plt.plot(conc_spot, conc_spot_im, 'ko', conc_spot, yth,'-k')
        if intercept>0:
            plt.text(conc_spot[nspot-1], conc_spot_im[1], 'y= %-3.4fx + %-3.4f\n' %(m,intercept))
        else:
            plt.text(conc_spot[nspot-1], conc_spot_im[1], 'y= %-3.4fx %-3.4f\n' %(m,intercept))
        plt.xlabel('pmol/mm2')
        plt.ylabel('drug/IS signal')
        plt.title('calibration curve')
        png1 = BytesIO()
        fig.savefig(png1, format='png')
        png2 = Image.open(png1)
        png2.save(sname + 'calibration_curve.tiff')
        plt.close()
        png1.close()
        mask_calibration_conc=(mask_drug-intercept)/(m*100)*mw_drug#concentrazione espressa come pg/pixel
        mask_calibration_conc=mask_calibration_conc*mask_tissue
        scale_factors=np.max(mask_calibration_conc)
        fig=plt.figure(dpi=300)
        plt.title(sname)
        if color_map:
            plt.imshow(mask_calibration_conc,aspect='equal',interpolation='None',cmap=color_map, vmin=0, vmax=scale_factors,origin='upper')
        else:
            plt.imshow(mask_calibration_conc,aspect='equal',interpolation='None', vmin=0, vmax=scale_factors,origin='upper')
        cbar=plt.colorbar()
        cbar.ax.set_xlabel(unit_F)
        plt.axis('off')
        png1 = BytesIO()
        fig.savefig(png1, format='png')
        png2 = Image.open(png1)
        png2.save(sname + '.tiff')
        plt.close()
        png1.close()
        coefficients=[(['m:']+[m]),(['q:']+[intercept]),(['weight:']+[weight_str])]
        calibration=[(['PTX(pmol)']+[conc_drug]),(['pixel_spot']+[dim_spot]),(['PTX/dPTX']+[conc_spot_im]),(['PTX(pmol/mm2)']+[conc_spot])]#corretto 
        with open('calibration.csv','wb')as csvfile:
            writer=csv.writer(csvfile)
            writer.writerows(calibration)
    return m,intercept,coefficients

def smooth_medianFilter(sname,MSImatrix_drug,MSImatrix_std,MSImatrix_tissue, mfiltrad,lob,scale_factors,r_p,c_p):
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
    val_piastra=np.percentile(nu, 99)   
    mask_tissue_bin= MSImatrix_tissue > val_piastra
    mask_tissue_bin=mask_tissue_bin.astype(int)

    np.savetxt(sname + 'mask_tissue_bin.msk',mask_tissue_bin,fmt ='%-3.2f',delimiter=',')    
    mask_drug = MSImatrix_drug * mask_tissue_bin
    mask_drug = ndimage.median_filter(mask_drug,mfiltrad) 
    mask_tissue = mask_drug> lob
    mask_tissue = mask_tissue.astype(int)
    mask_drug=mask_drug*mask_tissue # sogliata per il valore dei CTRL o di Otzu. 
    np.savetxt(sname + '.sim',mask_drug,fmt='%-3.2f',delimiter=',')
    np.savetxt(sname + '_drug.msk',mask_tissue,fmt ='%-3.2f',delimiter=',')
    return(mask_drug,mask_tissue,mask_tissue_bin)

def quantify(sname,mask_drug,mask_tissue,mask_tissue_bin,m,intercept,lod,mw_drug,voxel,th_s,scale_factors,unit_M, unit_F):   
    pixel_tot_tissue=float64(sum(mask_tissue_bin)) 
    pixel_neri=float64(sum(mask_tissue*mask_tissue_bin)) #come SM3
    percentuale=pixel_neri/pixel_tot_tissue*100
    mm2ROI=pixel_tot_tissue*voxel
    intensita_tot=float64(sum(mask_drug*mask_tissue_bin))
    drug_IS=intensita_tot/pixel_tot_tissue
    pmol=(drug_IS-intercept)/m*mm2ROI
    pmol_mg=pmol/(mm2ROI*th_s)
    conc=pmol_mg*mw_drug/1000
    mask_drug_conc=(mask_drug-intercept)/(m*100)*mw_drug#concentrazione espressa come pg/pixel
    mask_drug_conc=mask_drug_conc*mask_tissue_bin
    np.savetxt(sname + '_drug_conc.sim',mask_drug_conc,fmt ='%-3.2f',delimiter=',')
    lod_conc=(lod-intercept)/m
    return(pixel_tot_tissue,pixel_neri,percentuale,conc,mask_drug_conc,lod_conc)
    
def process(lob=[],lod=[],scale_factors=1,coord_spot=[],conc_drug=[],unit_M='ug/g',unit_F='ug/pixel',weight_type=[],mw_drug=[],voxel=0.01,th_s=0.01,r_p=4,c_p=4,color_map=[]):
    all_pixel_tot_tissue=['pixel_tot_tissue']
    all_pixel_neri=['Black_pixels']
    all_percentuale=['%']
    all_intensita_tot=[unit_M]
    vett_name=[]
    mfiltrad=3
    dircontent = os.listdir('.')
    selionimg = grep(dircontent,'.jpg')
    patterns_ctrl=['CTRL','ctrl','Ctrl']
    patterns_calibration=['CALIBRATION','calibration','Calibration']
    n_ctrl=0
    n_calibration=0 
    pos_calibration=[]
    for f in selionimg:
        a_r=np.zeros(3,dtype=bool)
        a_c=np.zeros(3,dtype=bool)
        i_c=0
        i_r=0
        sname = f[:-4] 
        for pattern_calibration in patterns_calibration:                      
            a_r[i_r]=re.search(pattern_calibration, sname)
            i_r=i_r+1
        if np.any(a_r):            
            pos_calibration=n_calibration
        n_calibration=n_calibration+1
        for pattern_ctrl in patterns_ctrl:           
            a_c[i_c]=re.search(pattern_ctrl, sname)
            i_c=i_c+1
        if np.any(a_c):
            pos_ctrl=n_ctrl
        n_ctrl=n_ctrl+1    
    if pos_ctrl > pos_calibration :
       if pos_ctrl-pos_calibration==1:
          if pos_calibration==0: 
              selionimg_ord=[selionimg[pos_calibration]]+[selionimg[pos_ctrl]]+selionimg[pos_ctrl+1:]
          else:
              selionimg_ord=[selionimg[pos_calibration]]+[selionimg[pos_ctrl]]+selionimg[:pos_calibration]+selionimg[pos_ctrl+1:]
       elif pos_ctrl==len(selionimg)-1:
           if pos_calibration ==0: 
               selionimg_ord=[selionimg[pos_calibration]]+[selionimg[pos_ctrl]]+selionimg[pos_calibration+1:pos_ctrl]
           else:
               selionimg_ord=[selionimg[pos_calibration]]+[selionimg[pos_ctrl]]+selionimg[:pos_calibration]+selionimg[pos_calibration+1:pos_ctrl]
       else:
           if pos_ctrl==0: 
               selionimg_ord=[selionimg[pos_calibration]]+[selionimg[pos_ctrl]]+selionimg[pos_calibration+1:pos_ctrl]+selionimg[pos_ctrl+1:]
           else:
               selionimg_ord=[selionimg[pos_calibration]]+[selionimg[pos_ctrl]]+selionimg[:pos_calibration]+selionimg[pos_calibration+1:pos_ctrl]+selionimg[pos_ctrl+1:]
    else:
       if pos_calibration-pos_ctrl==1:
          if pos_ctrl==0: 
              selionimg_ord=[selionimg[pos_calibration]]+[selionimg[pos_ctrl]]+selionimg[pos_calibration+1:]
          else:
              selionimg_ord=[selionimg[pos_calibration]]+[selionimg[pos_ctrl]]+selionimg[:pos_ctrl]+selionimg[pos_calibration+1:]
       elif pos_calibration==len(selionimg)-1:
           if pos_ctrl==0: 
               selionimg_ord=[selionimg[pos_calibration]]+[selionimg[pos_ctrl]]+selionimg[pos_ctrl+1:pos_calibration]
           else:
               selionimg_ord=[selionimg[pos_calibration]]+[selionimg[pos_ctrl]]+selionimg[:pos_ctrl]+selionimg[pos_ctrl+1:pos_calibration]
       else:
           if pos_ctrl==0: 
               selionimg_ord=[selionimg[pos_calibration]]+[selionimg[pos_ctrl]]+selionimg[pos_ctrl+1:pos_calibration]+selionimg[pos_calibration+1:]
           else:
               selionimg_ord=[selionimg[pos_calibration]]+[selionimg[pos_ctrl]]+selionimg[:pos_ctrl]+selionimg[pos_ctrl+1:pos_calibration]+selionimg[pos_calibration+1:]  
    print selionimg_ord 
    
    for f_ord in selionimg_ord:
        sname = f_ord[:-4]  
        print('PROCESSING %s : ' % (sname))
        MSImatrix_drug = np.genfromtxt(sname+'_MSImatrix_drug.sim',dtype=float,delimiter=',')
        MSImatrix_tissue=np.genfromtxt(sname+'_MSImatrix_tissue.sim',dtype=float,delimiter=',')
        MSImatrix_std=np.genfromtxt(sname+'_MSImatrix_std.sim',dtype=float,delimiter=',')
        mask_drug,mask_tissue,mask_tissue_bin=smooth_medianFilter(sname,MSImatrix_drug,MSImatrix_std,MSImatrix_tissue, mfiltrad,lob,scale_factors,r_p,c_p)
        a_r=np.zeros(3,dtype=bool)
        i_r=0
        for pattern_calibration in patterns_calibration:            
            a_r[i_r]=re.search(pattern_calibration, sname)
            i_r=i_r+1
        if np.any(a_r):
            m,intercept,coefficients=calibration_curve(sname,mask_drug,mask_tissue,coord_spot,conc_drug,weight_type,mw_drug,scale_factors,color_map,unit_F)
        
        else:
            pixel_tot_tissue,pixel_neri,percentuale,conc,mask_drug_conc,lod_conc=quantify(sname,mask_drug,mask_tissue,mask_tissue_bin,m,intercept,lod,mw_drug,voxel,th_s,scale_factors,unit_M, unit_F)
            all_pixel_tot_tissue=all_pixel_tot_tissue+[pixel_tot_tissue]
            all_pixel_neri=all_pixel_neri+[pixel_neri]
            all_percentuale=all_percentuale+[percentuale]
            all_intensita_tot=all_intensita_tot+[conc]
            figures_production(sname,mask_drug_conc,mask_tissue,mask_tissue_bin,scale_factors,color_map,voxel,lod_conc,unit_F)
        vett_name=vett_name+[sname]
    
    results=[(vett_name),(all_pixel_tot_tissue),(all_pixel_neri), (all_percentuale),(all_intensita_tot)]
    with open('results.csv','wb')as csvfile:
        writer=csv.writer(csvfile)
        writer.writerows(results)
        writer.writerows(coefficients)

def main():    
    parser = argparse.ArgumentParser(description="script to quantify drug concentration into tissue slice")
    parser.add_argument('-lob',dest = "lob",type = float,default=[],help = "LOB value used as threshold to obtain the drug_mask")  
    parser.add_argument('-lod',dest = "lod",type = float,default=[],help = "LOD value")  
    parser.add_argument('-sp',dest="coord_spot",type=int, nargs ='+', default=[], help=": localization of spots of the calibration curve into the tissue, i.e. the row and column numbers of the corners (four numbers ) for each spot")
    parser.add_argument('-c',dest="conc_drug",type=float, nargs ='+', default=[], help="drug concentration of calibration spots")
    parser.add_argument('-w',dest = "weight_type",type = float,default=[],help = "type of weight to fit calibration curve i.e 1: 1/y2 2: 1/x2 default empty no weight.")     
    parser.add_argument('-fs',dest = "scale_factors",type = float,default=1.0,help = "Maximum intensity in the color map. (Default fs=1.0)")
    parser.add_argument('-uf',dest="unit_F",type=str,  default='pg/pixel', help="unit of measure of the drug concentration in the figure, default value is pg/pixel")
    parser.add_argument('-um',dest="unit_M",type=str, default='ug/g', help="unit of measure of the mean drug concentration, default value is ug/g")
    parser.add_argument('-mw',dest="mw_drug",type=float, default=[], help="drug molecular weight i.e PTX 853")
    parser.add_argument('-v',dest="voxel",type=float, default=0.01, help=" pixel area (mm2)")
    parser.add_argument('-th_s',dest="th_s",type=float, default=0.01, help="tissue slice Thickness(mm)")
    parser.add_argument('-rp',dest="r_p",type=float, default=4, help="size (number of rows) of slice's edge used to estimate the plate intensity of the tissue ion peak")
    parser.add_argument('-cp',dest="c_p",type=float, default=4, help="size (number of columns) of slice's edge used to estimate the plate intensity of the tissue ion peak")
    parser.add_argument('-co',dest = "color_map",type = str,default=[],help = "color of ion.")     
    args = parser.parse_args()
    process(lob=args.lob,lod=args.lod,weight_type=args.weight_type,scale_factors=args.scale_factors,coord_spot=args.coord_spot,conc_drug=args.conc_drug,unit_F=args.unit_F,unit_M=args.unit_M,mw_drug=args.mw_drug,voxel=args.voxel,th_s=args.th_s,r_p=args.r_p,c_p=args.c_p,color_map=args.color_map)

if __name__ == '__main__':
   main()    