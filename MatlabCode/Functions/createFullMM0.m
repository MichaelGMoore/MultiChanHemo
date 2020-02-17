function Snew = createFullMM0(S)
% 10/102019
% Michael G. Moore, Michigan State University

% Create a meta-model from all sessions using the cross-validated lambda
% for L2-norm regularization

% This is intended for use in demixing a GCaMP mouse from the same line as
% a cohort of GFP mice

%%

% compute the meta-model coeffs
Ptrain = [];
Rtrain = [];
for s = 1:S.numSessions
    cond = S.session(s).PWR.FVE >= S.params.FVEThresh;
    Ptrain = cat(1,Ptrain,S.session(s).MM.predictors(cond,:));
    Rtrain = cat(1,Rtrain,S.session(s).MM.responses(cond,:));  
end
Rmean = mean(Rtrain,1);
Rtrain = reshape(Rtrain-Rmean,[size(Rtrain,1),S.numR*S.numF]);
Ptrain = zscore(Ptrain,1,1); % standardize the predictors

numP = size(Ptrain,2);
MMC = ((Ptrain'*Ptrain))\(Ptrain'*Rtrain);

S.CVMM0.numPredictors = numP;
S.CVMM0.ResponseMean = Rmean;
S.CVMM0.MMcoef = MMC;

Snew = S;

end

