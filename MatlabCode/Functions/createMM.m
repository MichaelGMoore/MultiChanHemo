function Snew = createMetaModels(S)
% 10/102019
% Michael G. Moore, Michigan State University


%% ********************************************************************
% CollectPredictors for Meta-Model from each session
% *********************************************************************
% Train the meta-model only on pixels with PWR-FVE > FVEThresh

% Use only original 19-predictors unless S.params.MM.morePredictors = true
  

for s = 1:S.numSessions

    N = S.session(s).numSamples;
    D = S.session(s).numPixels;
    numScales = S.params.MM.numScales;
   
    % collect the meta-model predictors
    if ~S.params.MM.morePredictors
        numP = 6*S.numR + .5*S.numR*(S.numR-1)  + numScales*S.numF;
    elseif S.params.MM.morePredictors
        numP = 6*S.numR + .5*S.numR*(S.numR-1)  + numScales*S.numF + 2*S.numR;
    end
        
    predictors = zeros(D,numP);
    p0 = 0;
    % L1 norm
    predictors(:,p0+1:S.numR) = S.session(s).stats.L1;
    p0 = p0+S.numR;
    % L1 squared
    predictors(:,p0+(1:S.numR)) = S.session(s).stats.L1.^2;
    p0 = p0+S.numR;
    % L2 
    predictors(:,p0+(1:S.numR)) = S.session(s).stats.P2.^.5;
    p0 = p0+S.numR;
    % L2^2
    predictors(:,p0+(1:S.numR)) = S.session(s).stats.P2.^1;
    p0 = p0+S.numR;
    % skew
    predictors(:,p0+(1:S.numR)) = S.session(s).stats.P3./(S.session(s).stats.P2.^1.5);
    p0 = p0+S.numR;
    % Kurtosis
    predictors(:,p0+(1:S.numR)) = S.session(s).stats.P4./(S.session(s).stats.P2.^2);
    p0 = p0+S.numR;
    % correlations between R channels
    for r1 = 1:S.numR
        for r2 = (r1+1):S.numR
            predictors(:,p0+r1+(r2-2)*S.numR) = squeeze(S.session(s).stats.PTP(r1,:,r2));
        end
    end
    p0 = p0 + .5*S.numR*(S.numR-1);
    
    % Blood Vessel Maps
    for f = 1:S.numF
        meanImg = S.session(s).meanFImg{f};
        for j = 1:numScales
            scale = 2^(j-1);
            meanImgBlurred = imgaussfilt(meanImg,scale);
            bvmap = meanImgBlurred./meanImg-1;           
            predictors(:,p0 + j) = bvmap(S.session(s).mask);
        end
        p0 = p0 + numScales;
    end
    
    % additional Blood Vessel Maps
    if S.params.MM.morePredictors
        for r = 1:S.numR
            meanImg = S.session(s).meanRImg{r};
            for j = 1:2
                scale = 2^j;
                meanImgBlurred = imgaussfilt(meanImg,scale);
                bvmap = meanImgBlurred./meanImg-1;           
                predictors(:,p0 +j) = bvmap(S.session(s).mask);
            end   
            p0 = p0 + 2;
        end
    end
    % save predictors
    S.session(s).MM.numP = numP;
    S.session(s).MM.predictors = predictors;
    S.session(s).MM.responses = S.session(s).PWR.coef;
end

% ************************************************************************
% compute the per-session meta-model coeffs and FVE

for s = 1:S.numSessions
    D = S.session(s).numPixels;
    cond = S.session(s).PWR.FVE >= S.params.FVEThresh;
    
    % train on the conditioned pixels
    R = S.session(s).MM.responses(cond,:);
    Rmean = mean(R,1);
    R = reshape(R-Rmean,[size(R,1),S.numR*S.numF]);
    P = S.session(s).MM.predictors(cond,:);
    P = zscore(P,1,1);
    C = P\R;
    
    % compute demixing coefficients for all pixels
    P = S.session(s).MM.predictors;
    P = zscore(P,1,1);    
    coef = reshape(Rmean + P*C,[D,S.numR,S.numF]); % [D x numR x numF];
   
    varI = zeros(D,S.numF);
    varF = zeros(D,S.numF);
    FVE = zeros(D,S.numF); 
    for d = 1:D
        PTP = permute(S.session(s).stats.PTP(:,d,:),[1,3,2]);
        PTR = permute(S.session(s).stats.PTR(:,d,:),[1,3,2]); 
        RTR = permute(S.session(s).stats.RTR(:,d,:),[1,3,2]);
       
        C = permute(coef(d,:,:),[2,3,1]); 
        varI(d,:) = diag(RTR);
        varF(d,:) = diag(RTR-2*C'*PTR+C'*PTP*C);
        FVE(d,:) = 1-varF(d,:)./varI(d,:); 
    end
    
    % save result
    S.session(s).MM.coef = coef;
    S.session(s).MM.varI = varI;
    S.session(s).MM.varF = varF;
    S.session(s).MM.FVE = FVE;   
end


Snew = S;
end

