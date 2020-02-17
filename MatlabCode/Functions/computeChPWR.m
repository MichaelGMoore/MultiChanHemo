function Snew = computeChPWR(S)
% 10/28/2019
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
    varI = zeros(D,S.numR,S.numF); % initial variance on each pixel
    varF = zeros(D,S.numR,S.numF); 
    FVE = zeros(D,S.numR,S.numF); % fractional variance explained by Pixelwise Regression
    % loop over pixels
    for d = 1:D 
        % loop over reflectance channels
        for r = 1:S.numR
            PTP = permute(S.session(s).stats.PTP(r,d,r),[1,3,2]);
            PTR = permute(S.session(s).stats.PTR(r,d,:),[1,3,2]); 
            RTR = permute(S.session(s).stats.RTR(:,d,:),[1,3,2]);

            C = PTP\PTR; % [1 x numF]
            coef(d,r,:) = C;  % [D x 1 x numF]
            varI(d,r,:) = diag(RTR);
            varF(d,r,:) = diag(RTR-2*C'*PTR+C'*PTP*C);
            FVE(d,r,:) = 1-varF(d,r,:)./varI(d,r,:); 
        end
    end
    S.session(s).ChPWR.coef = coef;   % [D x numR x numF] 
    S.session(s).ChPWR.varI = varI;   % [D x numF]
    S.session(s).ChPWR.varF = varF;   % [D x numF]
    S.session(s).ChPWR.FVE = FVE;     % [D x numF]
   
end

Snew = S;


end

