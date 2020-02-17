function Snew = computeBeerLambertCoefs(S)
% The input variable S contains the location of the HDF5 data files
% From these files we can pull out the spectra and compute Beer-Lambert
% demixing coefficients.

% pull the Prahl Extinction coeffs
load(fullfile(S.mainFolder,'BLdata','PrahlExtinctionCoeffs.mat'));

% pull the Hillman pathlengths
load(fullfile(S.mainFolder,'BLdata','HillmanPathlengths.mat'));

% loop over data files
for s = 1:S.numSessions
    % pull out the spectra
    spectra = readSpectra(S.fileNames{s},S.numF,S.numR);
    
    % Interpolate extinction coeffs and Hillman pathlengths onto spectra wavelengths
    %   multiply by .1 to convert from /cm to /mm
    E(:,1) = .1*interp1(extinction_coeffs.data(:,1),extinction_coeffs.data(:,2),spectra.lambda,'pchip','extrap');
    E(:,2) = .1*interp1(extinction_coeffs.data(:,1),extinction_coeffs.data(:,3),spectra.lambda,'pchip','extrap');
    X = interp1(hillman_pathlengths.data(:,1),hillman_pathlengths.data(:,2),spectra.lambda,'pchip','extrap');
    
    % compute the "characteristic wavelength approximation"
    for mu = 1:size(spectra.data,2)
        lambda(mu) = spectra.meanLambda(mu);
        x(mu) = interp1(spectra.lambda,X,lambda(mu),'makima','extrap');
        for nu = 1:2
            e(mu,nu) = interp1(spectra.lambda,E(:,nu),lambda(mu),'makima','extrap'); 
            M(mu,nu) = -x(mu)*e(mu,nu);     
        end
    end
    % we need to compress the excitation and emission (add .5*M_Fex + .5*M_Fem = M_F) 
    AF = repelem(.5*eye(S.numF),1,2);
    AR = eye(S.numR);
    M = blkdiag(AF,AR)*M; % dimensions of M are now (S.numF+S.numR) x 2
    % separate fluor and refl blocks
    MF = M(1:S.numF,:);
    MR = M((S.numF+1):end,:);
    C = MF/MR; % [numF x numR]
    C = C'; % transpose to get [numR x numF]

    % loop over pixels
    D = S.session(s).numPixels;
    coef = zeros(D,S.numR,S.numF);
    varI = zeros(D,S.numF);
    varF = zeros(D,S.numF);
    FVE = zeros(D,S.numF); 
    for d = 1:D 
        PTP = permute(S.session(s).stats.PTP(:,d,:),[1,3,2]);
        PTR = permute(S.session(s).stats.PTR(:,d,:),[1,3,2]); 
        RTR = permute(S.session(s).stats.RTR(:,d,:),[1,3,2]);

        coef(d,:,:) = C;  % [D x numR x numF]
        
        varI(d,:) = diag(RTR);
        varF(d,:) = diag(RTR-2*C'*PTR+C'*PTP*C);
        FVE(d,:) = 1-varF(d,:)./varI(d,:); 
    end
    S.session(s).BL0.coef = coef;   % [D x numR x numF] 
    S.session(s).BL0.varI = varI;   % [D x numF]
    S.session(s).BL0.varF = varF;   % [D x numF]
    S.session(s).BL0.FVE = FVE;
    
    
    % compute the "spectral Beer-Lambert"
    % choose reasonable estimates for background hemoglobing concentrations
    c(1,1) = (2.2e-3)*.85*.04; % oxygenated
    c(2,1) = (2.2e-3)*.15*.04; % deoxygenated
    w = [.5,.5,1,1]; % weights for path lengths
    for mu = 1:size(spectra.data,2)
        f = spectra.data(:,mu).*exp(-(w(mu)*X).*(E*c));
        for nu = 1:2
            M(mu,nu) = sum(f.*(-w(mu)*X).*E(:,nu))/sum(f);
        end
    
    end
    % we need to compress the excitation and emission (add .5*M_Fex + .5*M_Fem = M_F) 
    AF = repelem(eye(S.numF),1,2);
    AR = eye(S.numR);
    M = blkdiag(AF,AR)*M; % dimensions of M are now (S.numF+S.numR) x 2
    % separate fluor and refl blocks
    MF = M(1:S.numF,:);
    MR = M((S.numF+1):end,:);
    C = MF/MR; % [numF x numR]
    % compute demixing matrix (D0 is the crude approximation)
    C = C';   % transpose to get [numR x numF] 
    
    % loop over pixels
    D = S.session(s).numPixels;
    coef = zeros(D,S.numR,S.numF);
    varI = zeros(D,S.numF);
    varF = zeros(D,S.numF);
    FVE = zeros(D,S.numF); 
    for d = 1:D 
        PTP = permute(S.session(s).stats.PTP(:,d,:),[1,3,2]);
        PTR = permute(S.session(s).stats.PTR(:,d,:),[1,3,2]); 
        RTR = permute(S.session(s).stats.RTR(:,d,:),[1,3,2]);

        coef(d,:,:) = C;  % [D x numR x numF]
        
        varI(d,:) = diag(RTR);
        varF(d,:) = diag(RTR-2*C'*PTR+C'*PTP*C);
        FVE(d,:) = 1-varF(d,:)./varI(d,:); 
    end
    S.session(s).BL.coef = coef;   % [D x numR x numF] 
    S.session(s).BL.varI = varI;   % [D x numF]
    S.session(s).BL.varF = varF;   % [D x numF]
    S.session(s).BL.FVE = FVE;
    
end
clear s




Snew = S;
end

