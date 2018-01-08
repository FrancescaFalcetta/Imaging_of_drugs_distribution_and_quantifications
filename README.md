# Preprocessing_and_quantification_MSI
Imaging of drugs distribution and quantifications

This repository accompanies the article: Mass Spectrometry Imaging of drugs distribution in different tumour models: preprocessing and quantification issues.  It provides a reference for implementation of our pipeline for production of drug distribution and drug quantification images. At the moment, our pipeline accepts data from the single analytical platform, i.e. analyze 7.5.

There are four main functions in our pipeline:
a) getionimage  
b) intensity image production  
c) LOB and LOD calculation 
d) quantification
 
Figure 1.The proposed MSI data analysis pipeline
Reading raw data and preprocessing steps are implemented in the getionimage function.  The intensity image production function gives the drug intensity image and tissue mask. The quantification function plots the calibration curve,  quantifies the drug into the tissue and produces quantitative images.

Requirments
Python version 2.7
Python packages: numpy, sys, scipy, PIL, matplotlib, os, struct,skimage, argparse,
import itertools, pandas, re,rpy2, csv 
R version above 3.0
R packages: MALDIquantForeign, radiomics


Getionimage
Getionimage extracts ions  to be detected for each pixel as described in article. 

Inputs
-locIS: m/z range where to find the internal standard (IS) peak; 
-d_drug:  m/z difference between IS and drug peaks;
-locTiss: m/z range where to find the tissue peak;	
-hws and -hws_t:  half -width size of IS, drug peaks (-hws) and tissue peak (-hws_t, used only when diefferent from -hws);
-fr: Radius of the median filter
-th_mask: threshold value for the drug mask, given as input or automatically calculated from control sample. 
-fs: maximum intensity in the color map; 
-a:  modality to calculate the peak area. Options:   
	 a=[]: No background subtraction (method A in the article)
	-a range: trapezoidal background subtraction (method B in the article)
	-a min:   rectangular background subtraction (method C in the article)
-rp and -cp : size (number of rows and columns) of slice's edge used to estimate the plate intensity of the tissue ion peak; 

Example
python getionimage.py -locIS 289.06 289.34 -locTiss 281.20 281.50  -hws 3 -d_drug 5.03138 -fs 1.5 -a range

Outputs
Intensity MATRIX of drug, IS and tissue and a figure of drug intensity, tissue mask and drug mask.  
This figure is only to verify the goodness of ion extraction. Best figure will be generated with intensity image production.

intensity image production  
This script produces the images of drug intensity distribution, drug mask and tissue mask.

Inputs
-co: choice of color map type; 
-fs : maximum intensity in the color map 
-rp and -cp : size (number of rows and columns) of slice's edge used to estimate the plate intensity of the tissue ion peak;
-th_mask: threshold value for the drug mask, given as input or automatically calculated from control sample. 

Example
python intensity_image_production.py -co inferno -fs 1.5

Outputs
Intensity MATRIX of drug, IS and tissue and a figure of drug intensity, figure of drug mask, and figure of tissue mask.

   
lob and lod calculation 
This script calculates the LOD and LOB in a calibration curve sample.
 As described in the article: The LOB is defined as the signal at the 95th percentile of replicate blank measures. The LOD is the mean signal of a sample drug concentration where only the 5% of the observed replicated values is under the LOB. 

The localization of spots can be found opening the calibration_drug.msk file and/or the calibation.sim file with excel® or another program to identify the start and end rows and columns as shown in the figure.

 

Inputs
-sp: localization of spots of the calibration curve into the tissue, i.e. the row and column numbers of the corners (four numbers ) for each spot;  
-c : drug concentration of calibration spots;
The order of localization and concentration of spot must be the same.

Example
python lob_and_lod_calculation.py -sp 6 21 27 43 25 40 28 43 50 65 23 41 11 25 7 22 30 43 4 18 -c 10 5 2.5 1 0.5 0
The file name of calibration curve must contain “calibration” inside it (not sensitive case).

Outputs
LOB e LOD value printed on video.


Quantification image production
Script to quantify the drug concentration into the samples. 

The localization of spots can be found opening the calibration_drug.msk file and/or the calibation.sim file with excel® or another program to identify the start and end rows and columns as shown in the figure as show in the example below.

Inputs
-lob: LOB value used as threshold to obtain the drug_mask;
-lod: LOD value;
-sp: localization of spots of the calibration curve into the tissue, i.e. the row and column numbers of the corners (four numbers ) for each spot;  
-c : drug concentration of calibration spots;
The order of localization and concentration of spot must be the same.
-w: type of weight to fit calibration curve:  
-w 1: 1/y2
-w 2: 1/x2  
-co: choice of color map type; 
-fs : maximum intensity in the color map 
-rp and -cp : size (number of rows and columns) of slice's edge used to estimate the plate intensity of the tissue ion peak;
-mw: drug molecular weight (i.e 853 for PTX);
-um: unit of measure of the mean drug concentration, default value is ug/g;
-uf: unit of measure of the drug concentration in the figure, default value is pg/pixel;
-v: pixel area (mm2);
-th_s: tissue slice thickness (mm).

Example
python quantification.py -lob 0.0120 -lod 0.24  -w 1 -mw 853 -fs 5 -co inferno -sp 6 21 27 43 25 40 28 43 50 65 23 41 11 25 7 22 30 43 4 18 -c 10 5 2.5 1 0.5 0 
The file name of calibration curve must contain “calibration” inside it (not sensitive case).

Output 
Drug quantitative images and a table reporting the number of pixels of the tissue, of the number of pixels above the Lob in the drug image and the mean drug concentration.
 

Contact
email: Francesca Falcetta (francesca.falcetta@marionegri.it)

