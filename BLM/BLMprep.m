% Script to prepare data elements for Beer Lambert Models:
% 20180614

%% ************************************************************************
%  BEGIN BEER-LAMBERT MODEL PREPARATION SECTION
%  ************************************************************************

%% Model parameters
%   http://rstb.royalsocietypublishing.org/content/royptb/suppl/2016/08/23/rstb.2015.0360.DC1/rstb20150360supp1.pdf
% or from correspondence with Matt Valley

% We have pathlength estimates from 400-700nm in 2nm increments
% We have extinction coeffs from 250-1000 nm in 2nm increments

% our spectral distributions run from 380-780 nm in 1 nm increments. 

% To cover the path-lengths from 700-780 nm we will have to extrapolate the tabulated data

opsBLM=struct; % Beer-Lambert Model params structure

% make a wavelengths table running from 400-800 nm in 2 nm increments
opsBLM.lambdas = (380:780)';

% extract Scott Prahl's Hemoglobin extinction coeffs
load PrahlExtinctionCoeffs.mat
k = find(extinction_coeffs.data(:,1) >= 380 & extinction_coeffs.data(:,1) <= 780);
opsBLM.E(:,1:2) = interp1(extinction_coeffs.data(k,1),extinction_coeffs.data(k,2:3),opsBLM.lambdas,'linear','extrap');
clear k extinction_coeffs
% convert from per-centimeter to per-millimeter
opsBLM.E = .1*opsBLM.E;

% extract corresponding pathlengths
load Hilman_pathlengths.mat
% extrapolate to find remaining values (extrapolated results look very reasonable)
opsBLM.x = interp1(hilman_pathlengths.data(:,1),hilman_pathlengths.data(:,2),opsBLM.lambdas,'pchip','extrap');
clear k hilman_pathlengths

% load spectral distributions
spectra = load('MattValleySpectra.mat');
k = find(spectra.data(:,1) >= 380 & spectra.data(:,1) <= 780);

% extract relevant spectra from Matt Valley's Measurements:
opsBLM.PS1 = spectra.data(k,4);
opsBLM.PS1 = opsBLM.PS1/sum(opsBLM.PS1);
opsBLM.PS2 = spectra.data(k,5);
opsBLM.PS2 = opsBLM.PS2/sum(opsBLM.PS2);
opsBLM.PSex = spectra.data(k,2);
opsBLM.PSex = opsBLM.PSex/sum(opsBLM.PSex);

% include Matt Valleys recorded emission spectra
% Notes:  ,'linewidth',2
%   These have passed throught tissue and are presumably modified by absorption
%   We have an alternative source from the Tseih lab
opsBLM.PFemMV1C = spectra.data(k,3);
opsBLM.PFemMV1C = opsBLM.PFemMV1C/sum(opsBLM.PFemMV1C);
opsBLM.PFemMV2C = spectra.data(k,6);
opsBLM.PFemMV2C = opsBLM.PFemMV2C/sum(opsBLM.PFemMV2C);

clear k spectra

% load emission filter distribution for Semrock FF01-525/45 Optical Filter
fname = 'MattValleyEmissionFilter';
filterData=dlmread(fname,'\t',4,0);
opsBLM.PCem = interp1(filterData(:,1),filterData(:,2),opsBLM.lambdas,'pchip','extrap');
clear fname filterData

% for the Fluorphore spectra, we can use the Tsein data
fname = 'TseinFluorophore_Excitation_Spectra.xlsx';
FexData = xlsread(fname); % column 1 is wavelength and column 5 is GFP data
FexData(isnan(FexData))=0; % set missing values to 0
opsBLM.PFex = interp1(FexData(:,1),FexData(:,5),opsBLM.lambdas,'pchip','extrap');
opsBLM.PFex = opsBLM.PFex/sum(opsBLM.PFex);

fname = 'TseinFluorophore_Emission_Spectra.xlsx';
FemData = xlsread(fname); % column 1 is wavelength and column 5 is GFP data
FemData(isnan(FemData))=0; % set missing values to zero
opsBLM.PFem = interp1(FemData(:,1),FemData(:,5),opsBLM.lambdas,'pchip','extrap');
opsBLM.PFem = opsBLM.PFem/sum(opsBLM.PFem);

clear fname FexData FemData

%% save results

save('/home/mmoore/Dropbox/Work/projects/WFOM_BOLD_3color/BLM/opsBLM.mat','opsBLM')

%% check distributions
figure
hold on
plot(opsBLM.lambdas,opsBLM.PSex/max(opsBLM.PSex),'b','linewidth',2)
plot(opsBLM.lambdas,opsBLM.PS1/max(opsBLM.PS1),'y','linewidth',2)
plot(opsBLM.lambdas,opsBLM.PS2/max(opsBLM.PS2),'r','linewidth',2)
plot(opsBLM.lambdas,opsBLM.PFex/max(opsBLM.PFex),'c','linewidth',2)
plot(opsBLM.lambdas,opsBLM.PFem/max(opsBLM.PFem),'g','linewidth',2)
plot(opsBLM.lambdas,opsBLM.PFemMV1C/max(opsBLM.PFemMV1C),'m--','linewidth',2)
plot(opsBLM.lambdas,opsBLM.PFemMV2C/max(opsBLM.PFemMV2C),'m.-','linewidth',2)
plot(opsBLM.lambdas,opsBLM.PCem/max(opsBLM.PCem),'k','linewidth',2)

%% check pathlengths
figure
plot(opsBLM.lambdas,opsBLM.x,'k')

%% check estinction coeffs
figure
plot(opsBLM.lambdas,opsBLM.E(:,1),'r',opsBLM.lambdas,opsBLM.E(:,2),'b')
 
%% Super-pose pathlengths and distributions

figure
hold on
plot(opsBLM.lambdas,opsBLM.PSex/max(opsBLM.PSex),'b','linewidth',2)
plot(opsBLM.lambdas,opsBLM.PFem/max(opsBLM.PFem),'g','linewidth',2)
plot(opsBLM.lambdas,opsBLM.PS1/max(opsBLM.PS1),'y','linewidth',2)
plot(opsBLM.lambdas,opsBLM.PS2/max(opsBLM.PS2),'r','linewidth',2)
plot(opsBLM.lambdas,opsBLM.x/max(opsBLM.x),'k','linewidth',2)

legend('F_{ex}','F_{em}','B_1','B_2','pathlength')

%% super-pose distributions onto absorbance factors
figure
hold on
area(opsBLM.lambdas,opsBLM.PSex.*opsBLM.PFex/max(opsBLM.PSex.*opsBLM.PFex),'FaceColor','b')
area(opsBLM.lambdas,opsBLM.PCem.*opsBLM.PFem/max(opsBLM.PCem.*opsBLM.PFem),'FaceColor','g')
area(opsBLM.lambdas,opsBLM.PS1/max(opsBLM.PS1),'FaceColor',[1 .7 0])
area(opsBLM.lambdas,opsBLM.PS2/max(opsBLM.PS2),'FaceColor','r')
plot(opsBLM.lambdas,opsBLM.E(:,1).*opsBLM.x/max(opsBLM.E(:,2).*opsBLM.x),'k-','Linewidth',2)
plot(opsBLM.lambdas,opsBLM.E(:,2).*opsBLM.x/max(opsBLM.E(:,2).*opsBLM.x),'k--','Linewidth',2)
xlabel('wavelength[nm]')
yticks([])
legend('Excitation','Emission','575 Reflectiance','630 Reflectance','HbO Absorbance','HbR Absorbance')
title('Hemoglobin Absorbence and Optical Spectra')
box on
set(gcf,'Color','w')