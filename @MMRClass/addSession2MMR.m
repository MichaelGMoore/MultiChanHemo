function obj = addSession2MMR(obj,data)
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
    clear m
    
    obj.sessionID{mSession} = data.sessionID;
    obj.animalID{mSession} = data.animalID;
    obj.folder{mSession} = data.folder;
    obj.fname{mSession} = data.fname;
        
    obj.mask{mSession} = data.mask;

    y = data.y;
    D = sum(data.mask(:));
    N = data.frames;
    
    % parse the reflectance and fluorescence
    reflCh = find(contains(data.typeCh,'R'));
    fluoCh = find(contains(data.typeCh,'F'));
    numReflCh = length(reflCh);
    numFluoCh = length(fluoCh);
    labelR = cell(numReflCh,1); 
    labelF = cell(1,numFluoCh); 
    for nr = 1:numReflCh
        obj.MMdata{mSession}.labelR{nr,1} = [num2str(data.waveCh(reflCh(nr))),'nm'];
    end
    for nf = 1:numFluoCh
        obj.MMdata{mSession}.labelF{1,nf} = [num2str(data.waveCh(fluoCh(nf))),'nm'];
    end

    % First compute the Meta Model Responses
    for d=1:D
        R = [];
        % build the predictor array (N x numReflCh)
        for r = reflCh
            R = cat(2,R,data.y{r}(:,d));
        end
        for nf = 1:numFluoCh
            f = fluoCh(nf);
            F = data.y{f}(:,d);    
            coef = R\F;
            varI = var(F,1,1);
            varF = var(F-R*coef,1,1);
            obj.MMdata{mSession}.FVE(d,nf) = 1-varF/varI; 
            obj.MMdata{mSession}.responses(d,:,nf) = coef;
        end  
    end
    clear d r nf R f F coef
        
    % compute meta model predictors
    obj.MMdata{mSession}.predictors_labels = {};
    obj.MMdata{mSession}.predictors = [];    
    % L-1 norm of each reflectance channel
    for nr = 1:numReflCh
        r = reflCh(nr);
        obj.MMdata{mSession}.predictors_labels{end+1,1} = [labelR{nr},' L1'];    
        temp = mean(abs(y{r}),1);
        obj.MMdata{mSession}.predictors = cat(1,obj.MMdata{mSession}.predictors,temp); 
    end   
    % L-1 norm squared of each reflectance channel
    for nr = 1:numReflCh
        r = reflCh(nr);
        obj.MMdata{mSession}.predictors_labels{end+1,1} = [labelR{nr},' L1^2'];    
        temp = mean(abs(y{r}),1);
        obj.MMdata{mSession}.predictors = cat(1,obj.MMdata{mSession}.predictors,temp.^2); 
    end    
    % L-2 norm of each reflectance channel
    for nr = 1:numReflCh
        r = reflCh(nr);
        obj.MMdata{mSession}.predictors_labels{end+1,1} = [labelR{nr},' L2'];    
        temp = mean(y{r}.^2,1);
        obj.MMdata{mSession}.predictors = cat(1,obj.MMdata{mSession}.predictors,temp.^.5); 
    end
     % L-2 norm squared of each reflectance channel
    for nr = 1:numReflCh
        r = reflCh(nr);
        obj.MMdata{mSession}.predictors_labels{end+1,1} = [labelR{nr},' L2^2'];    
        temp = mean(y{r}.^2,1);
        obj.MMdata{mSession}.predictors = cat(1,obj.MMdata{mSession}.predictors,temp); 
    end   
    % skewness of each reflectance channel
    for nr = 1:numReflCh
        r = reflCh(nr);
        obj.MMdata{mSession}.predictors_labels{end+1,1} = [labelR{nr},' skewness'];    
        temp1 = mean(y{r}.^3,1);
        temp2 = mean(y{r}.^2,1);
        obj.MMdata{mSession}.predictors = cat(1,obj.MMdata{mSession}.predictors,temp1./(temp2.^1.5)); 
    end    
    % kurtosis of each reflectance channel
    for nr = 1:numReflCh
        r = reflCh(nr);
        obj.MMdata{mSession}.predictors_labels{end+1,1} = [labelR{nr},' kurtosis'];    
        temp1 = mean(y{r}.^4,1);
        temp2 = mean(y{r}.^2,1);
        obj.MMdata{mSession}.predictors = cat(1,obj.MMdata{mSession}.predictors,temp1./(temp2.^2)); 
    end
    % correlations between reflectance channels
    for nr1 = 1:numReflCh
        r1 = reflCh(nr1);
        for nr2 = (nr1+1):numReflCh
            r2 = reflCh(nr2);
            obj.MMdata{mSession}.predictors_labels{end+1,1} = [labelR{nr1},'-',labelR{nr2},' correlation']; 
            temp = mean(y{r1}.*y{r2},1);
            obj.MMdata{mSession}.predictors = cat(1,obj.MMdata{mSession}.predictors,temp); 
        end
    end
    % blood vessel maps
    num_scales = obj.params.MMR.numScales;
    f = fluoCh(1);
    meanImg = data.meanImg{f};
    mask = data.mask;
    for j=1:num_scales
        obj.MMdata{mSession}.predictors_labels{end+1,1} = ['BV map ',num2str(2^(j-1))];
        scale = 2^(j-1);
        meanImgBlurred = imgaussfilt(meanImg,scale);
        temp = meanImgBlurred./meanImg-1;           
        temp = temp(mask)';
        obj.MMdata{mSession}.predictors = cat(1,obj.MMdata{mSession}.predictors,temp); 
    end      
    
    % if everything runs without errors, update numSession
    if mSession > obj.numSessions
        obj.numSessions = obj.numSessions + 1;       
    end
end

