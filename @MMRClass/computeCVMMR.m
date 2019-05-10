function obj = computeCVMMR(obj,data)
% compute demixing statistics for each animal using only training data from
% other animals in the same cre line

for s = 1:obj.numSessions
    test1 = isequal(obj.sessionID{s},data.sessionID);
    test2 = isequal(obj.animalID{s},data.animalID);
    if test1 && test2
        break
    end
end

% parse the reflectance and fluorescence
reflCh = find(contains(data.typeCh,'R'));
fluoCh = find(contains(data.typeCh,'F'));
numReflCh = length(reflCh);
numFluoCh = length(fluoCh);
labelR = cell(numReflCh,1); 
labelF = cell(1,numFluoCh); 
for nr = 1:numReflCh
    obj.CVMMR{s}.labelR{nr,1} = [num2str(data.waveCh(reflCh(nr))),'nm'];
end
for nf = 1:numFluoCh
    obj.CVMMR{s}.labelF{1,nf} = [num2str(data.waveCh(fluoCh(nf))),'nm'];
end

FVEThresh = obj.params.MMR.FVEThresh;
obj.CVMMR{s}.FVEThresh = FVEThresh;

y = data.y;
D = sum(data.mask(:));
N = data.frames;

% train the model using good pixels from other animals only
% we need responses and predictors for each reflectance channel

for nf = 1:numFluoCh
    Responses = [];
    Predictors = [];
    for s2 = 1:obj.numSessions
        if ~isequal(obj.animalID{s2},data.animalID)
            cond = obj.MMdata{s2}.FVE(:,nf) > FVEThresh;
            r = obj.MMdata{s2}.responses(cond,:,nf);
            Responses = cat(2,Responses,r');
            p = obj.MMdata{s2}.predictors(:,cond); 
            Predictors = cat(2,Predictors,p);         
        end   
    end
    RespMean = mean(Responses,2);
    Predictors = zscore(Predictors,1,2);
    Coeffs = (Responses - RespMean)/Predictors;
    
    obj.CVMMR{s}.RespMean(:,nf) = RespMean;
    obj.CVMMR{s}.Coeffs(:,:,nf) = Coeffs;  
    
    Predictors = obj.MMdata{s}.predictors; % for demixing we need to use all pixels
    coef = RespMean + Coeffs*zscore(Predictors,1,2);
    
    f = fluoCh(nf);
    for d = 1:D
        R = [];
        for r = reflCh
            R = cat(2,R,y{r}(:,d));
        end
        F = y{f}(:,d);
        varI = var(F,1,1);
        varF = var(F-R*coef(:,d),1,1);
        obj.CVMMR{s}.varI(d,nf) = varI;
        obj.CVMMR{s}.varF(d,nf) = varF;
        obj.CVMMR{s}.FVE(d,nf) = 1-varF/varI; 
    end
end




end

