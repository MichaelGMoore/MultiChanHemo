function Snew = createFullMM(S)
% 10/102019
% Michael G. Moore, Michigan State University

% Create a meta-model from all sessions using the cross-validated lambda
% for L2-norm regularization

% This is intended for use in demixing a GCaMP mouse from the same line as
% a cohort of GFP mice

%%
lambda = S.CVMM.lambdaOpt;

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
MMC = ((Ptrain'*Ptrain)+lambda^2*eye(numP))\(Ptrain'*Rtrain);

S.CVMM.numPredictors = numP;
S.CVMM.ResponseMean = Rmean;
S.CVMM.MMcoef = MMC;

Snew = S;

end

