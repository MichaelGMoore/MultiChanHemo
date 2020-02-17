function Snew = createCVMM(S)
% 10/102019
% Michael G. Moore, Michigan State University

% Create the leave-one-animal out cross-validated Meta-model with L2
% regularization

% find the unique set of animals
animalList = unique({S.session.animalID});
numAn = length(animalList);


lamList = exp(linspace(log(1e-2),log(1e3),70));
lamList(1) = 0;  

cost = zeros(size(lamList)); % this will be the final variance after demixing
for l = 1:length(lamList)
    lambda = lamList(l);
       
    % loop over animals
    costAn = zeros(numAn,1);
    for an = 1:numAn
        % train the model using only predictors from the other animals
        % train only on pixels satisfying PWR.FVE > FVEThresh
        Ptrain = [];
        Rtrain = [];
        for s = 1:S.numSessions
            if ~isequal(animalList{an},S.session(s).animalID)
                cond = S.session(s).PWR.FVE >= S.params.FVEThresh;
                Ptrain = cat(1,Ptrain,S.session(s).MM.predictors(cond,:));
                Rtrain = cat(1,Rtrain,S.session(s).MM.responses(cond,:));
            end
        end
        clear s
        Rmean = mean(Rtrain,1);
        Rtrain = reshape(Rtrain-Rmean,[size(Rtrain,1),S.numR*S.numF]);
        Ptrain = zscore(Ptrain,1,1); % standardize the predictors
        
        % compute the meta-model coefficient with L2 regularization
        numP = size(Ptrain,2);
        MMC = ((Ptrain'*Ptrain)+lambda^2*eye(numP))\(Ptrain'*Rtrain); % [numP x numR*numF] this is the meta-model coefficient
        
        % for each session of same animal, demix using MMC
        costAnPix = [];
        for s = 1:S.numSessions
            D = S.session(s).numPixels;
            if isequal(animalList{an},S.session(s).animalID)
                Ptest = S.session(s).MM.predictors;
                Ptest = zscore(Ptest,1,1); % standardize the predictors
                coef = reshape(Rmean + Ptest*MMC,[D,S.numR,S.numF]);
                % loop over pixels to compute FVE
                varF = zeros(D,S.numF);
                varI = zeros(D,S.numF);
                for d = 1:D
                    PTP = permute(S.session(s).stats.PTP(:,d,:),[1,3,2]);
                    PTR = permute(S.session(s).stats.PTR(:,d,:),[1,3,2]); 
                    RTR = permute(S.session(s).stats.RTR(:,d,:),[1,3,2]);

                    C = permute(coef(d,:,:),[2,3,1]); 
                    varI(d,:) = diag(RTR);
                    varF(d,:) = diag(RTR - 2*C'*PTR + C'*PTP*C);
                end
                clear d
                costAnPix = cat(1,costAnPix,varF); % [D x numF], per-pixel cost function
            end
            
        end
        clear s
        costAn(an) = median(costAnPix,[1,2]); % averaging over pixels treats each animal equally rather than bias towards an animal with more pixels
    end
    clear an
    % compute the total cost accross animals for the given lambda
    cost(l) = mean(costAn);
  
end
clear l

%% determine the optimal lambda

[~,lMin] = min(cost);
lambdaOpt = lamList(lMin);
S.CVMM.lambdaOpt = lambdaOpt;

%% Demix each session using meta-model without regularization

for s1 = 1:S.numSessions
    an1 = S.session(s1).animalID;
    % create training data from the other animals
    Ptrain = [];
    Rtrain = [];
    for s2 = 1:S.numSessions
        an2 = S.session(s2).animalID;
        if ~isequal(an1,an2)
            cond = S.session(s2).PWR.FVE >= S.params.FVEThresh;
            Ptrain = cat(1,Ptrain,S.session(s2).MM.predictors(cond,:));
            Rtrain = cat(1,Rtrain,S.session(s2).MM.responses(cond,:));        
        end
    end
    clear s2
    Rmean = mean(Rtrain,1);
    Rtrain = reshape(Rtrain-Rmean,[size(Rtrain,1),S.numR*S.numF]);
    Ptrain = zscore(Ptrain,1,1); % standardize the predictors
    numP = size(Ptrain,2);
    % compute the meta-model coefficient with no regularization
    MMC = (Ptrain'*Ptrain)\(Ptrain'*Rtrain);
    
    % demix using the umregularized MMC
    D = S.session(s1).numPixels;
        
    Ptest = S.session(s1).MM.predictors;           
    Ptest = zscore(Ptest,1,1); % standardize the predictors    
    coef = reshape(Rmean + Ptest*MMC,[D,S.numR,S.numF]); % form the coefficients
    % loop over pixels to compute FVE
    varI = zeros(D,S.numF);
    varF = zeros(D,S.numF);
    FVE = zeros(D,S.numF);
    for d = 1:D
        PTP = permute(S.session(s1).stats.PTP(:,d,:),[1,3,2]);
        PTR = permute(S.session(s1).stats.PTR(:,d,:),[1,3,2]); 
        RTR = permute(S.session(s1).stats.RTR(:,d,:),[1,3,2]);

        C = permute(coef(d,:,:),[2,3,1]); 
        varI(d,:) = diag(RTR);
        varF(d,:) = diag(RTR-2*C'*PTR+C'*PTP*C);
        FVE(d,:) = 1 - varF(d,:)./varI(d,:);
    end
    clear d
    % save results
    S.session(s1).CVMM0.MMC = MMC;
    S.session(s1).CVMM0.coef = coef;
    S.session(s1).CVMM0.varI = varI;
    S.session(s1).CVMM0.varF = varF;
    S.session(s1).CVMM0.FVE = FVE;   
end
clear s1

%% Demix each session using meta-model constructed from other animals with the optimal lambda

for s1 = 1:S.numSessions
    an1 = S.session(s1).animalID;
    % create training data from the other animals
    Ptrain = [];
    Rtrain = [];
    for s2 = 1:S.numSessions
        an2 = S.session(s2).animalID;
        if ~isequal(an1,an2)
            cond = S.session(s2).PWR.FVE >= S.params.FVEThresh;
            Ptrain = cat(1,Ptrain,S.session(s2).MM.predictors(cond,:));
            Rtrain = cat(1,Rtrain,S.session(s2).MM.responses(cond,:));        
        end
    end
    clear s2
    Rmean = mean(Rtrain,1);
    Rtrain = reshape(Rtrain-Rmean,[size(Rtrain,1),S.numR*S.numF]);
    Ptrain = zscore(Ptrain,1,1); % standardize the predictors
    numP = size(Ptrain,2);
    % compute the meta-model coefficient with L2 regularization
    MMC = ((Ptrain'*Ptrain)+lambdaOpt^2*eye(numP))\(Ptrain'*Rtrain);
    
    % demix using the regularized MMC
    D = S.session(s1).numPixels;
        
    Ptest = S.session(s1).MM.predictors;           
    Ptest = zscore(Ptest,1,1); % standardize the predictors    
    coef = reshape(Rmean + Ptest*MMC,[D,S.numR,S.numF]); % form the coefficients
    % loop over pixels to compute FVE
    varI = zeros(D,S.numF);
    varF = zeros(D,S.numF);
    FVE = zeros(D,S.numF);
    for d = 1:D
        PTP = permute(S.session(s1).stats.PTP(:,d,:),[1,3,2]);
        PTR = permute(S.session(s1).stats.PTR(:,d,:),[1,3,2]); 
        RTR = permute(S.session(s1).stats.RTR(:,d,:),[1,3,2]);

        C = permute(coef(d,:,:),[2,3,1]); 
        varI(d,:) = diag(RTR);
        varF(d,:) = diag(RTR-2*C'*PTR+C'*PTP*C);
        FVE(d,:) = 1 - varF(d,:)./varI(d,:);
    end
    clear d
    % save results
    S.session(s1).CVMM.MMC = MMC;
    S.session(s1).CVMM.coef = coef;
    S.session(s1).CVMM.varI = varI;
    S.session(s1).CVMM.varF = varF;
    S.session(s1).CVMM.FVE = FVE;   
end
clear s1

Snew = S;
end

