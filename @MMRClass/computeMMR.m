function obj = computeMMR(obj,data)
% compute demixing statistics for each animal using a meta model trained
% only on that animal

% find the dataset

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
    obj.MMR{s}.labelR{nr,1} = [num2str(data.waveCh(reflCh(nr))),'nm'];
end
for nf = 1:numFluoCh
    obj.MMR{s}.labelF{1,nf} = [num2str(data.waveCh(fluoCh(nf))),'nm'];
end

% condition on FVE > FVEThresh
FVEThresh = obj.params.MMR.FVEThresh;
obj.MMR{s}.FVEThresh = FVEThresh;

y = data.y;
D = sum(data.mask(:));
N = data.frames;

% we want metamodel coeffs for each fluo channel, because "responses" will
% be different. Predictors, based on refl channels, are the same in each
% case, but "good" pixels for training data may differ
for nf = 1:numFluoCh
    Responses = [];
    Predictors = [];
    cond = obj.MMdata{s}.FVE(:,nf) > FVEThresh;
    r = obj.MMdata{s}.responses(cond,:,nf);
    Responses = cat(2,Responses,r');
    p = obj.MMdata{s}.predictors(:,cond); % for training the model we use only the "good" pixels
    Predictors = cat(2,Predictors,p);

    RespMean = mean(Responses,2);
    Predictors = zscore(Predictors,1,2);
    Coeffs = (Responses - RespMean)/Predictors;
    
    obj.MMR{s}.RespMean(:,nf) = RespMean;
    obj.MMR{s}.Coeffs(:,:,nf) = Coeffs;
        
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
        obj.MMR{s}.varI(d,nf) = varI;
        obj.MMR{s}.varF(d,nf) = varF;
        obj.MMR{s}.FVE(d,nf) = 1-varF/varI; 
    end
end





end

