function obj = trainMM(obj)

% parse the reflectance and fluorescence
obj.trainedMM.labelR = obj.MMdata{1}.labelR;
obj.trainedMM.labelF = obj.MMdata{1}.labelF;
numReflCh = size(obj.trainedMM.labelR,1);
numFluoCh = size(obj.trainedMM.labelF,1);


FVEThresh = obj.params.MMR.FVEThresh;
obj.trainedMM.FVEThresh = FVEThresh;

% sum over fluo channels
for nf = 1:numFluoCh
    Responses = [];
    Predictors = [];
    % include all sessions
    for s2 = 1:obj.numSessions
        cond = obj.MMdata{s2}.FVE(:,nf) > FVEThresh;
        r = obj.MMdata{s2}.responses(cond,:,nf);
        Responses = cat(2,Responses,r');
        p = obj.MMdata{s2}.predictors(:,cond); 
        Predictors = cat(2,Predictors,p);          
    end
    RespMean = mean(Responses,2);
    Predictors = zscore(Predictors,1,2);
    Coeffs = (Responses - RespMean)/Predictors;
    
    obj.trainedMM.RespMean(:,nf) = RespMean;
    obj.trainedMM.Coeffs(:,:,nf) = Coeffs;  
    
end



end

