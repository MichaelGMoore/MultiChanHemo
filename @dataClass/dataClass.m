classdef dataClass
    % class to store multichannel data for hemo demixing
    % each instance of the class is a single recording session
    
    properties

        folder      % location of raw data
        fname       % name of raw data file
        
        animalID        % identifier of animal
        sessionID   % identifier of recording session       
        numCh       % number of channels
        nameCh      % name of channel
        waveCh      % approximate wavelength of each channel in [nm]
        typeCh      % 'F' or 'R' for fluorescence/reflectance
        
        fs          % sampling frequency
        frames      % array of frame indices
        times       % array of time stamps
        
        meanImg     % mean images of each channel
        mask        % mask 
        y           % delta I over mean I for each channel (channels ordered by wavelength)
        
        filtered    % if empty, not filtered, otherwise a frequency and width indicated
        timeShifts  % if empty, not shifted, otherwise a vector of frame-shifts applied to each channel
                
    end
    
    methods
        function obj = dataClass(images,mask,metaData)
            
            tempSize = size(images{1});
            assert(length(tempSize)==3,'images must be 3D array')
            for k = 2:length(images)
                assert(isequal(size(images{k}),tempSize),'all channel data must have same dimensions')
            end
            clear k tempSize
            
            if isfield(metaData,'folder')
                obj.folder = metaData.folder;
            end
            if isfield(metaData,'fname')
                obj.fname = metaData.fname;
            end
            
            if isfield(metaData,'animalID')
                obj.animalID = metaData.animalID;
            end
            if isfield(metaData,'sessionID')
                obj.sessionID = metaData.sessionID;
            end
            
            if isfield(metaData,'numCh') 
                assert(metaData.numCh == length(images),'unexpected number of channels found')
                obj.numCh = metaData.numCh;
            else
                obj.numCh = length(images);
            end
            if isfield(metaData,'nameCh')
                obj.nameCh = metaData.nameCh;
            end
            if isfield(metaData,'waveCh')
                obj.waveCh = metaData.waveCh;
            end
            if isfield(metaData,'typeCh')
                obj.typeCh = metaData.typeCh;
            else
                obj.typeCh{1} = 'F';
                for k = 2:length(images)
                    obj.typeCh{k} = 'R';
                end
                disp('assuming first channel is fluorescence')
            end
            
            if isfield(metaData,'fs')
                obj.fs = metaData.fs;
            else
                error('must provide sampling frequency')
            end
            obj.frames = size(images{1},1);
            obj.times = [1:obj.frames]'/obj.fs;
            
            for k = 1:obj.numCh
                obj.meanImg{k} = squeeze(mean(images{k},1));
            end
            obj.mask = mask;
            for k = 1:obj.numCh
                num = double(images{k}(:,mask));
                denom = obj.meanImg{k}(mask)';
                obj.y{k} = num./denom -1;
            end         
            
        end
        
        function ch = getFluoChan(obj)
            ch= find(contains(obj.typeCh,'F'));
    
        end
        
        % external methods
        inspectData(obj)
        inspectTimeShifts(obj,timeShifts)
        inspectLinearity(obj)
        obj = applyTimeShifts(obj,params)
        obj = applyLowPassFilter(obj,params);
       
        
    end
end

