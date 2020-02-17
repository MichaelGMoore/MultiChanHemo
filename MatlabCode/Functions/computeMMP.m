function predictors = computeMMP(stats,meanFImg,meanRImg,mask,S)

% 10/30/2019
% Michael G. Moore, Michigan State University

% compute the predictors to apply the meta-model


D = sum(mask(:));
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
predictors(:,p0+1:S.numR) = stats.L1;
p0 = p0+S.numR;
% L1 squared
predictors(:,p0+(1:S.numR)) = stats.L1.^2;
p0 = p0+S.numR;
% L2 
predictors(:,p0+(1:S.numR)) = stats.P2.^.5;
p0 = p0+S.numR;
% L2^2
predictors(:,p0+(1:S.numR)) = stats.P2.^1;
p0 = p0+S.numR;
% skew
predictors(:,p0+(1:S.numR)) = stats.P3./(stats.P2.^1.5);
p0 = p0+S.numR;
% Kurtosis
predictors(:,p0+(1:S.numR)) = stats.P4./(stats.P2.^2);
p0 = p0+S.numR;
% correlations between R channels
for r1 = 1:S.numR
    for r2 = (r1+1):S.numR
        predictors(:,p0+r1+(r2-2)*S.numR) = squeeze(stats.PTP(r1,:,r2));
    end
end
p0 = p0 + .5*S.numR*(S.numR-1);

% Blood Vessel Maps
for f = 1:S.numF
    meanImg = meanFImg{f};
    for j = 1:numScales
        scale = 2^(j-1);
        meanImgBlurred = imgaussfilt(meanImg,scale);
        bvmap = meanImgBlurred./meanImg-1;           
        predictors(:,p0 + j) = bvmap(mask);
    end
    p0 = p0 + numScales;
end

% additional Blood Vessel Maps
if S.params.MM.morePredictors
    for r = 1:S.numR
        meanImg = meanRImg{r};
        for j = 1:2
            scale = 2^j;
            meanImgBlurred = imgaussfilt(meanImg,scale);
            bvmap = meanImgBlurred./meanImg-1;           
            predictors(:,p0 +j) = bvmap(mask);
        end   
        p0 = p0 + 2;
    end
end


end

