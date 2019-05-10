% script to load and process all 3 cre-lines of GFP data through the
% following Regression Models:

% Pixel-Wise Regression
    % linReg        all reflectance channels
    % linReg1Ch     single reflectance channels
    
% Beer-Lambert Regression
    % BLR           single-wavelength Beer-Lambert model
    % SBLR          integrated-wavelength Beer-Lambert model
    
% Meta Model regression
    % MMR           training and testing meta-model on single animal
    % CVMMR         training data on other animals and testing it on leave-one-out animal
    % trainedMM     model trained on all animals in creline, for use on GCaMP data
       
creNames = {'Cux2','Ntsr1','Rorb'};

for c = 1:length(creNames)
    creName = creNames{c};
    disp(creName)
    % create a list of the Cux2 GFP dataset filenames and locations
    folder = ['/home/mmoore/Optical_Data_Analysis/Data/Widefield/Matt_Valley/GFP_Mice/',creName,'/Raw'];
    fdir = dir(folder);
    fdir = fdir(3:end);  % first 2 entries are . and ..

    fileList = struct;
    for f = 1:length(fdir)
        fileList(f,1).folder = folder;
        fileList(f,1).fname = fdir(f).name;
    end
    clear f fdir folder

    % Pre-Processing Parameters

    params = struct;

    % file data
    params.fileList = fileList;

    % Low Pass Filter Parameters
    params.filter.f0 = 5;  % filter frequency in Hz
    params.filter.width = 0; % width of the filter edge in Hz

    % time shift parameters
    params.timeShifts = zeros(length(fileList),3);
    if isequal(creName,'Cux2')
        params.timeShifts(3,:) = [0,4,2];
    end

    % PWR params
    params.PWR = [];

    % BLR params
    params.BLR.BLRcoefs = [1.031361189458454;0.128940249873746];  % basic Beer-Lambert coefficients
    params.BLR.SBLRcoefs = [1.176023650535313;-0.434043363539705]; % spectral Beer-Lambert coefficients

    % MMR params
    params.MMR.numScales = 6; % number of scales for blood vessel maps
    params.MMR.FVEThresh = .8; % threshold on PWR FVE for inclusing in MM training data

    % create  various regression objects for the cohort
    PWR = PWRClass(creName,params); % Pixel-Wise Regression
    BLR = BLRClass(creName,params); % Beer-Lambert Regression
    MMR = MMRClass(creName,params); % Meta-Model Regression

    % Load the data sets as needed, using a custom parser to convert to dataClass object
    % select a file
    for f = 1:length(fileList)
        disp(['adding session ',num2str(f)])

        % call a custom parser to read each file
        folder = fileList(f).folder;
        fname = fileList(f).fname;
        clear images mask metaData
        [images,mask,metaData] = readMattValleyh5(folder,fname);
        clear folder fname

        % Use images, mask, and metaData to create a dataClass object
        clear data
        data = dataClass(images,mask,metaData);
        clear images mask metaData

        % Pre-process the data:
            % Filter data (optional)
            data = applyLowPassFilter(data,params);
            % Re-align data (optional)
            data = applyTimeShifts(data,params);

        % add session to the PixelWiseRegression Class
        PWR = addSession2PWR(PWR,data);

        % add session to the BeerLambertRegression Class
        BLR = addSession2BLR(BLR,data);   

        % add session to the MetaModeRegression Class
        MMR = addSession2MMR(MMR,data);

    end

    % Run the meta-model analysis
    % these have to be run after all sessions are added to the metamodel class
    for f = 1:length(fileList)
        disp(['analyzing Meta Models for ',num2str(f)])
        % call a custom parser to read each file
        folder = fileList(f).folder;
        fname = fileList(f).fname;
        clear images mask metaData
        [images,mask,metaData] = readMattValleyh5(folder,fname);
        clear folder fname

        % Use images, mask, and metaData to create a dataClass object
        clear data
        data = dataClass(images,mask,metaData);
        clear images mask metaData

        % Pre-process the data:
            % Filter data (optional)
            data = applyLowPassFilter(data,params);
            % Re-align data (optional)
            data = applyTimeShifts(data,params);

        % evaluate demixing efficiency of meta-model
        MMR = computeMMR(MMR,data);

        % cross-validation test of meta-model
        MMR = computeCVMMR(MMR,data);

    end

    % Train the full Cre-line MetaModel

    MMR =  trainMM(MMR);

    % save results 
    save(['SavedObjects/PWR',creName],'PWR')
    save(['SavedObjects/BLR',creName],'BLR')
    save(['SavedObjects/MMR',creName],'MMR')

end



































