% Pre-processing datasets for MultiChannel Hemodynamics demixing

% 10/16/2019
% Michael G. Moore, Michigan State University

% Run this script first to read raw data and extract necessary components
% for training and evaluating the various demixing models


%% specify project parameters
% edit this section to reflect your dataset

% this will be the data-structure containing results for the cohort
S = struct; 

% specify the path to the code
S.mainFolder = '/home/mmoore/Dropbox/Work/projects/WFOM_MultiChanHemo/Version3';
addpath(genpath(S.mainFolder));

% specify number of fluoresence and reflectance channels:
S.numF = 1;    % number of fluoresence channels
S.numR = 2;     % number of reflectance channels 

% give channels names
S.chName = {'520 nm','577 nm','630 nm'};

% specify the sampling frequency of the data in Hertz
S.Fs = 100;

% Specify low-pass filter parameters (if desired):
S.params.LPF.doLPF = true;    % set true to low-pass-filter the data
S.params.LPF.f0 = 5;   % filter frequency in Hz

S.params.FVEThresh = .7; % threshold on PWR FVE for inclusion in study

% MMR params
S.params.MM.numScales = 6; % number of scales for blood vessel maps
S.params.MM.morePredictors = true;

%% Pre-process the Cux2 dataset

% % assign a name to the dataset
% S.dataSetName = 'Cux2';   
% % give path to .h5 files of imaging sessions:
% S.dataFolder = '/home/mmoore/Data/MultichanHemoData/Cux2'; 
% 
% analyzePixels

%% Pre-process the Rorb dataset

% assign a name to the dataset
S.dataSetName = 'Rorb';   
% give path to .h5 files of imaging sessions:
S.dataFolder = '/home/mmoore/Data/MultichanHemoData/Rorb'; 

analyzePixels

%% Pre-process the Ntsr1 dataset

% assign a name to the dataset
S.dataSetName = 'Ntsr1';   
% give path to .h5 files of imaging sessions:
S.dataFolder = '/home/mmoore/Data/MultichanHemoData/Ntsr1'; 

analyzePixels
