
%% Run the pipeline

% Steps:

%   2   read the spectra and compute Beer Lambert coefficients
%   3   compute the spatial model and cross-validation

%% specify path to the project
addpath(genpath('/home/mmoore/Dropbox/Work/projects/WFOM_MultiChanHemo/Version3'));

%% choose a dataset 
S = SNtsr1;

%% run the pixelwise regression
S = computePWR(S);

%% run single-channel regressions
S = computeChPWR(S);

%% run single-channel ratiometric
S = computeRatiometric(S);

%% use the spectral data to compute beer lambert coefficients
S = computeBeerLambertCoefs(S);

%% compute the base meta-model:
S = createMM(S);

%% compute the cross-valided meta-model
% should optimize for FVE or varF?
% we can choose regularization based on min-max or min-mean?

S = createCVMM(S);

%% Create the full meta-model for the dataset
S = createFullMM(S);

%% Create the full meta-model for the dataset without regularization
S = createFullMM0(S);

%% Rename the output and write to file

outName = ['S',S.dataSetName];
assignin('base',outName,S);
save(outName,outName)





















