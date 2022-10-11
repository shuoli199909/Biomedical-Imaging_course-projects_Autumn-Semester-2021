% Task 1.2(PET)
clear; close all; clc;

% Check "InvivoData.mat"
invivo_loc = "InvivoData.mat";
load(invivo_loc);
matrix = invivo.matrix;
figure(1);
DisplayData(invivo.petsino,[2,2,1]); title("PET-sinogram");
DisplayData(invivo.ctsino,[2,2,2]); title("CT-sinogram");

% Reconstruct both the PET and CT data
filter = CalcFilter(matrix+1);
invivo.petfbp = CalcFBPRecon(invivo.angles,matrix,invivo.petsino,filter);
invivo.ctfbp = CalcFBPRecon(invivo.angles,matrix,invivo.ctsino,filter);
DisplayData(invivo.petfbp,[2,2,3]); title("PET-FBP");
DisplayData(invivo.ctfbp,[2,2,4]); title("CT-FBP");

% Calculate linear attenuation coefficients from HU data
mue_water = 1.707E-01;                     % water  @ 100 keV  [/cm]
invivo.mue = invivo.ctfbp*mue_water/1000;