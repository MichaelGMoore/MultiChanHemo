function obj = addSession2PWR(obj,data)
    % adds a session to a pixelWiseRegression object
    
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
    
    % compute the pixelwise regression coefficients
    D = sum(data.mask(:));
    N = data.frames;

    % list the reflectance and fluorescence channels
    reflCh = find(contains(data.typeCh,'R'));
    fluoCh = find(contains(data.typeCh,'F'));

    numReflCh = length(reflCh);
    numFluoCh = length(fluoCh);

    % compute linear regression coefficients
    obj.linReg{mSession} = struct;
    obj.linReg{mSession}.coef = zeros(D,length(reflCh),length(fluoCh));

    % attach the wavelengths of the channels as labels
    obj.linReg{mSession}.labelR = cell(length(reflCh),1); 
    obj.linReg{mSession}.labelF = cell(1,length(fluoCh)); 
    for nr = 1:numReflCh
        obj.linReg{mSession}.labelR{nr,1} = [num2str(data.waveCh(reflCh(nr))),'nm'];
    end
    for nf = 1:numFluoCh
        obj.linReg{mSession}.labelF{1,nf} = [num2str(data.waveCh(fluoCh(nf))),'nm'];
    end

    % regression done pixelwise as advertised
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
            coef = R\F;
            obj.linReg{mSession}.coef(d,:,nf) = coef;
            obj.linReg{mSession}.varI(d,nf) = var(F,1,1);
            obj.linReg{mSession}.varF(d,nf) = var(F-R*coef,1,1);
            obj.linReg{mSession}.FVE(d,nf) = 1-obj.linReg{mSession}.varF(d,nf)/obj.linReg{mSession}.varI(d,nf);
        end  
    end
    
    % compute linear regression 1channel coefficients
    obj.linReg1Ch{mSession} = struct;
    obj.linReg1Ch{mSession}.coef = zeros(D,length(reflCh),length(fluoCh));
    % attach the wavelengths of the channels as labels
    obj.linReg1Ch{mSession}.labelR = cell(length(reflCh),1); 
    obj.linReg1Ch{mSession}.labelF = cell(1,length(fluoCh)); 
    for nr = 1:numReflCh
        obj.linReg1Ch{mSession}.labelR{nr,1} = [num2str(data.waveCh(reflCh(nr))),'nm'];
    end
    for nf = 1:numFluoCh
        obj.linReg1Ch{mSession}.labelF{1,nf} = [num2str(data.waveCh(fluoCh(nf))),'nm'];
    end

    % regression done pixelwise as advertised
    for d=1:D
        for nr = 1:numReflCh
            r = reflCh(nr);
            R = data.y{r}(:,d);
            % compute coeffs for each fluor channel
            for nf = 1:numFluoCh
                f = fluoCh(nf);
                F = data.y{f}(:,d);    
                coef = R\F;
                obj.linReg1Ch{mSession}.coef(d,nr,nf) = coef;
                obj.linReg1Ch{mSession}.varI(d,nr,nf) = var(F,1,1);
                obj.linReg1Ch{mSession}.varF(d,nr,nf) = var(F-R*coef,1,1);
                obj.linReg1Ch{mSession}.FVE(d,nr,nf) = 1 -obj.linReg1Ch{mSession}.varF(d,nr,nf)/obj.linReg1Ch{mSession}.varI(d,nr,nf);
            end  
        end
    end
       
    % if everything runs without errors, update numSession
    if mSession > obj.numSessions
        obj.numSessions = obj.numSessions + 1;       
    end
end

