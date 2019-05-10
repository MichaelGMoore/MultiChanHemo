%% Create the output file

outfile = ['H5Objects',filesep,'fullAnalysis'];

delete(outfile);

creNames = {'Cux2','Ntsr1','Rorb'};

for c = 1:3
    creName = creNames{c};
    load(['SavedObjects',filesep,'PWR',creName])
    load(['SavedObjects',filesep,'BLR',creName])
    load(['SavedObjects',filesep,'MMR',creName])
    
    for s = 1:PWR.numSessions
        datasetName = ['/',creName,'/Session',num2str(s)];
        % session info 
        sessionID = str2num(PWR.sessionID{s});
        animalID = str2num(PWR.animalID{s}(2:end));    
        h5put(outfile,datasetName,sessionID)
        h5put(outfile,datasetName,animalID)
        % pre-processing info
        filterF0_Hz = PWR.params.filter.f0;
        filterWidth_Hz = PWR.params.filter.width;
        frameShifts = PWR.params.timeShifts(s,:);  
        mask = double(PWR.mask{s});
        h5put(outfile,datasetName,filterF0_Hz)
        h5put(outfile,datasetName,filterWidth_Hz)
        h5put(outfile,datasetName,frameShifts)
        h5put(outfile,datasetName,mask)
        % initial variance
        varI = PWR.linReg{s}.varI;
        h5put(outfile,datasetName,varI)
        % PWR info        
        linRegcoef = PWR.linReg{s}.coef;
        linRegvarF = PWR.linReg{s}.varF;
        linRegFVE = PWR.linReg{s}.FVE;
        linReg1Chcoef = PWR.linReg1Ch{s}.coef;
        linReg1ChvarF = PWR.linReg1Ch{s}.varF;
        linReg1ChFVE = PWR.linReg1Ch{s}.FVE;
        h5put(outfile,datasetName,linRegcoef)
        h5put(outfile,datasetName,linRegvarF)
        h5put(outfile,datasetName,linRegFVE)
        h5put(outfile,datasetName,linReg1Chcoef)
        h5put(outfile,datasetName,linReg1ChvarF)
        h5put(outfile,datasetName,linReg1ChFVE)
        % BLR info
        BLRcoef = BLR.BLR{s}.coef;
        BLRvarF = BLR.BLR{s}.varF;
        BLRFVE = BLR.BLR{s}.FVE;
        SBLRcoef = BLR.SBLR{s}.coef;
        SBLRvarF = BLR.SBLR{s}.varF;
        SBLRFVE = BLR.SBLR{s}.FVE;
        h5put(outfile,datasetName,BLRcoef)
        h5put(outfile,datasetName,BLRvarF)
        h5put(outfile,datasetName,BLRFVE)
        h5put(outfile,datasetName,SBLRcoef)
        h5put(outfile,datasetName,SBLRvarF)
        h5put(outfile,datasetName,SBLRFVE)       
        % MMR info
        MMpredictors = MMR.MMdata{s}.predictors;
        MMpredictorsLabels = MMR.MMdata{s}.predictors_labels;
        MMRresponseMean = MMR.MMR{s}.RespMean;
        MMRcoefs = MMR.MMR{s}.Coeffs;
        MMRvarF = MMR.MMR{s}.varF;
        MMRFVE = MMR.MMR{s}.FVE;
        CVMMRresponseMean = MMR.CVMMR{s}.RespMean;
        CVMMRcoefs = MMR.CVMMR{s}.Coeffs;
        CVMMRvarF = MMR.CVMMR{s}.varF;
        CVMMRFVE = MMR.CVMMR{s}.FVE;
        h5put(outfile,datasetName,MMpredictors)
%         h5put(outfile,datasetName,MMpredictorsLabels)
        h5put(outfile,datasetName,MMRresponseMean)
        h5put(outfile,datasetName,MMRcoefs)
        h5put(outfile,datasetName,MMRvarF)
        h5put(outfile,datasetName,MMRFVE)
        h5put(outfile,datasetName,CVMMRresponseMean)
        h5put(outfile,datasetName,CVMMRcoefs)
        h5put(outfile,datasetName,CVMMRvarF)
        h5put(outfile,datasetName,CVMMRFVE)
    end
    % trained metamodel
    datasetName = ['/',creName,'/trainedMM'];
    trainedMMrespMean = MMR.trainedMM.RespMean;
    trainedMMcoef = MMR.trainedMM.Coeffs;
    h5put(outfile,datasetName,trainedMMrespMean)
    h5put(outfile,datasetName,trainedMMcoef)
end





%% *************************************************************************

% mask = double(data.mask);
% sessionID = data.sessionID;
% animalID = data.animalID;
% filtered = data.filtered;
% timeShifts = data.timeShifts;
% 
% outfile = ['./MidlineDemixedData/',sessionID,'_',animalID,'_MM'];
% 
% h5create(outfile,'/Data/filtered',size(filtered),'ChunkSize',[1,size(filtered,2)],'DataType','double');
% h5write(outfile,'/Data/filtered',filtered)
% 
% h5create(outfile,'/Data/timeShifts',size(timeShifts),'ChunkSize',[1,size(timeShifts,2)],'DataType','double');
% h5write(outfile,'/Data/timeShifts',timeShifts)
% 
% h5create(outfile,'/Data/mask',size(mask),'ChunkSize',[1,size(mask,2)],'DataType','double');
% h5write(outfile,'/Data/mask',mask)
% 
% h5create(outfile,'/Data/F_PWR',size(F_PWR),'ChunkSize',[1,size(F_PWR,2)],'DataType','double');
% h5write(outfile,'/Data/F_PWR',F_PWR)
% 
% h5create(outfile,'/Data/F_MM2',size(F_MM2),'ChunkSize',[1,size(F_MM2,2)],'DataType','double');
% h5write(outfile,'/Data/F_MM2',F_MM2)
% 
% h5create(outfile,'/Data/F_MM12',size(F_MM12),'ChunkSize',[1,size(F_MM12,2)],'DataType','double');
% h5write(outfile,'/Data/F_MM12',F_MM12)