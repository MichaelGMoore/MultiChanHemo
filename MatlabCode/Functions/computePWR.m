function Snew = computePWR(S)
% 10/09/2019
% Michael G. Moore, Michigan State University

% *************************************************************************
% use the pixel-wise statistics to compute the pixel-wise regression
% *************************************************************************

% this computes the PWR coef and FVE  without regularization

%% loop over sessions
for s = 1:S.numSessions
    
    N = S.session(s).numSamples;
    D = S.session(s).numPixels;

    coef = zeros(D,S.numR,S.numF); % optimal PWR coeffs for each pixel
    varI = zeros(D,S.numF); % initial variance on each pixel
    varF = zeros(D,S.numF); 
    FVE = zeros(D,S.numF); % fractional variance explained by Pixelwise Regression
    % loop over pixels
    for d = 1:D 
        PTP = permute(S.session(s).stats.PTP(:,d,:),[1,3,2]);
        PTR = permute(S.session(s).stats.PTR(:,d,:),[1,3,2]); 
        RTR = permute(S.session(s).stats.RTR(:,d,:),[1,3,2]);

        C = PTP\PTR; % [numR x numF]
        coef(d,:,:) = C;  % [D x numR x numF]
        
        varI(d,:) = diag(RTR);
        varF(d,:) = diag(RTR-2*C'*PTR+C'*PTP*C);
        FVE(d,:) = 1-varF(d,:)./varI(d,:); 
    end
    S.session(s).PWR.coef = coef;   % [D x numR x numF] 
    S.session(s).PWR.varI = varI;   % [D x numF]
    S.session(s).PWR.varF = varF;   % [D x numF]
    S.session(s).PWR.FVE = FVE;
   
end

Snew = S;


end

