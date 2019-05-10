function obj = addSession2BLR(obj,data)
    % adds a session to a BeerLambertRegression object
    
    % each field is a cell-array with an entry for each session
    
    % data is a dataClass object for the session
    assert(isa(data,'dataClass'),'data must be a dataClass object')

    % Assume session not included, set mSession = number sessions + 1
    mSession = obj.numSessions + 1;
    % check if session already exists, and if it does, set mSession =
    % session number to overwrite it.   
    for m = 1:length(obj.sessionID)
        sesID = obj.sessionID{m};
        anID = obj.animalID{m};
        if isequal(data.sessionID,sesID) && isequal(data.animalID,anID)
            warning('This Dataset already exists in this pixelWiseRegression object. The analysis will be overwritten.')
            mSession = m;
        end
    end
    
    obj.sessionID{mSession} = data.sessionID;
    obj.animalID{mSession} = data.animalID;
    obj.folder{mSession} = data.folder;
    obj.fname{mSession} = data.fname;
        
    obj.mask{mSession} = data.mask;

    D = sum(data.mask(:));
    N = data.frames;

    % list the reflectance and fluorescence channels
    reflCh = find(contains(data.typeCh,'R'));
    fluoCh = find(contains(data.typeCh,'F'));

    numReflCh = length(reflCh);
    numFluoCh = length(fluoCh);

    obj.BLR{mSession} = struct;
    coef = obj.params.BLR.BLRcoefs;

    % attach the wavelengths of the channels as labels
    obj.BLR{mSession}.labelR = cell(length(reflCh),1); 
    obj.BLR{mSession}.labelF = cell(1,length(fluoCh)); 
    for nr = 1:numReflCh
        obj.BLR{mSession}.labelR{nr,1} = [num2str(data.waveCh(reflCh(nr))),'nm'];
    end
    for nf = 1:numFluoCh
        obj.BLR{mSession}.labelF{1,nf} = [num2str(data.waveCh(fluoCh(nf))),'nm'];
    end
    
    % regress with spatially uniform map
    for d=1:D
        R = [];
        % build the predictor array (N x numReflCh)
        for r = reflCh
            R = cat(2,R,data.y{r}(:,d));
        end
        % compute coeffs for each fluor channel
        for nf = 1:numFluoCh
            f = fluoCh(nf);
            F = data.y{f}(:,d);    
            obj.BLR{mSession}.coef(d,:,nf) = coef(:,nf);
            obj.BLR{mSession}.varI(d,nf) = var(F,1,1);
            obj.BLR{mSession}.varF(d,nf) = var(F-R*coef(:,nf),1,1);
            obj.BLR{mSession}.FVE(d,nf) = 1-obj.BLR{mSession}.varF(d,nf)/obj.BLR{mSession}.varI(d,nf);
        end  
    end
    
    obj.SBLR{mSession} = struct;
    coef = obj.params.BLR.SBLRcoefs;

    % attach the wavelengths of the channels as labels
    obj.SBLR{mSession}.labelR = cell(length(reflCh),1); 
    obj.SBLR{mSession}.labelF = cell(1,length(fluoCh)); 
    for nr = 1:numReflCh
        obj.SBLR{mSession}.labelR{nr,1} = [num2str(data.waveCh(reflCh(nr))),'nm'];
    end
    for nf = 1:numFluoCh
        obj.SBLR{mSession}.labelF{1,nf} = [num2str(data.waveCh(fluoCh(nf))),'nm'];
    end
    
    % regress with spatially uniform map
    for d=1:D
        R = [];
        % build the predictor array (N x numReflCh)
        for r = reflCh
            R = cat(2,R,data.y{r}(:,d));
        end
        % compute coeffs for each fluor channel
        for nf = 1:numFluoCh
            f = fluoCh(nf);
            F = data.y{f}(:,d);    
            obj.SBLR{mSession}.coef(d,:,nf) = coef(:,nf);
            obj.SBLR{mSession}.varI(d,nf) = var(F,1,1);
            obj.SBLR{mSession}.varF(d,nf) = var(F-R*coef(:,nf),1,1);
            obj.SBLR{mSession}.FVE(d,nf) = 1-obj.SBLR{mSession}.varF(d,nf)/obj.SBLR{mSession}.varI(d,nf);
        end  
    end
    
    % if everything runs without errors, update numSession
    if mSession > obj.numSessions
        obj.numSessions = obj.numSessions + 1;       
    end
end

