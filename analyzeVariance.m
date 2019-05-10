% set the path to the code base
addpath(genpath('/home/mmoore/Dropbox/Work/projects/WFOM_MultiChanHemo/Version2'))

%% load the data

creName = {'Cux2','Ntsr1','Rorb'};

clear medianBLR medianMMR medianCVMMR
for c = 1:length(creName)
    load(['BLR' creName{c}]);
    load(['MMR' creName{c}]);
    S = BLR.numSessions; % number of sessions
    % get a list of animal IDs
    clear anID
    for s = 1:S
        anID{s} = BLR.animalID{s};
    end
%     clear s
    [uA,ia,ic]  = unique(anID,'stable');
    A = length(uA); % number of animals
    % combine data for unique animal IDs
    for a = 1:A
        varI = [];
        BLRvarF = [];
        MMRvarF = [];
        CVMMRvarF = [];
        for s = 1:S
            if ic(s) == a
                varI = cat(1,varI,BLR.BLR{s}.varI);
                BLRvarF = cat(1,BLRvarF,BLR.SBLR{s}.varF);
                MMRvarF = cat(1,MMRvarF,MMR.MMR{s}.varF);
                CVMMRvarF = cat(1,CVMMRvarF,MMR.CVMMR{s}.varF);
            end
        end
%         clear s
        BLR_FVR = BLRvarF./varI;
        MMR_FVR = MMRvarF./varI;
        CVMMR_FVR = CVMMRvarF./varI;
       
        medianBLR(c,a) = median(BLR_FVR);
        medianMMR(c,a) = median(MMR_FVR);
        medianCVMMR(c,a) = median(CVMMR_FVR);
    end
%     clear a
end
% clear c
medianBLR(medianBLR==0) = NaN;
medianMMR(medianMMR==0) = NaN;
medianCVMMR(medianCVMMR==0) = NaN;

clearvars -except creName medianBLR medianMMR medianCVMMR

%%  take mean and std across animals
creNameC = creName';
SBL = nanmean(medianBLR,2);
MM = nanmean(medianMMR,2);
CVMM = nanmean(medianCVMMR,2);

mean = table(creNameC,SBL,MM,CVMM);

SBL = nanstd(medianBLR,1,2);
MM = nanstd(medianMMR,1,2);
CVMM = nanstd(medianCVMMR,1,2);

std = table(creNameC,SBL,MM,CVMM);


