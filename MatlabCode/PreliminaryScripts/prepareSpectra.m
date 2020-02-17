% Prepare the spectra needed to run Beer Lambert on Matt Valley GFP data.

spectra = load('MattValleySpectra.mat');
% column | data
% -------+-----------
%    1   | frequency   
%    2   | excitation LED
%    3   | emission (1 cam)
%    4   | reflectance 1
%    5   | reflectance 2
%    6   | emission (2 cam)

figure
hold on
for m = 2:6
    plot(spectra.data(:,1),spectra.data(:,m))
end
clear m
legend('ex','em(1 cam)','r1','r2','em(2 cam)')
title('Matt Valley Measured Spectra')

filterData=dlmread('MattValleyEmissionFilter','\t',4,0);
% this is the Semrock GFP filter data
figure
plot(filterData(:,1),filterData(:,2))
legend('Semrock GFP filter Spectra')

% for the Fluorphore spectra, we can use the Tsein data
FexData = xlsread('TseinFluorophore_Excitation_Spectra.xlsx'); % column 1 is wavelength and column 5 is GFP data
FexData(isnan(FexData))=0; % set missing values to 0
FemData = xlsread('TseinFluorophore_Emission_Spectra.xlsx'); % column 1 is wavelength and column 5 is GFP data
FemData(isnan(FemData))=0; % set missing values to zero
figure
plot(FexData(:,1),FexData(:,5),FemData(:,1),FemData(:,5))
legend('Ex','Em')
title('Tsein GFP spectra')

% put into schema
R1 = cat(2,spectra.data(:,1),spectra.data(:,4));
R2 = cat(2,spectra.data(:,1),spectra.data(:,5));
Fex = cat(2,spectra.data(:,1),spectra.data(:,2)...
    .*interp1(FexData(:,1),FexData(:,5),spectra.data(:,1),'pchip','extrap'));
Fem = cat(2,spectra.data(:,1),interp1(filterData(:,1),filterData(:,2),spectra.data(:,1),'pchip','extrap')...
    .*interp1(FemData(:,1),FemData(:,5),spectra.data(:,1),'pchip','extrap'));

figure
hold on
plot(R1(:,1),R1(:,2),'-','color',[1 .7 0],'linewidth',2)
plot(R2(:,1),R2(:,2),'-','color',[1 0 0],'linewidth',2)
plot(Fex(:,1),Fex(:,2),'-','color',[0 0 1],'linewidth',2)
plot(Fem(:,1),Fem(:,2),'-','color',[0 1 0],'linewidth',2)

clear FemData FexData filterData spectra

spectra = struct;
spectra.R1 = R1;
spectra.R2 = R2;
spectra.Fex = Fex;
spectra.Fem = Fem;

clear R1 R2 Fex Fem

save mySpectra.mat spectra

clear spectra