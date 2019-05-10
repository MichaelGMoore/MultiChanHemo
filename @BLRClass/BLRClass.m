classdef BLRClass
    %BLRCLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cohortName          % a name for the PWR cohort, e.g. 'Cux2'       

        params              % all pipeline parameters used in processing the data
        
        numSessions         % number of sessions included in the PWR cohort
                            %   updated dynamically as each session is added
        
        % cell array properties (cell labels each session)
        sessionID           % ID for each session
        animalID            % ID for each animal
        folder              % folder containing raw data
        fname               % raw data filename
        
        mask                % mask for each animal k
                                                       
        BLR               % regression using Beer-Lambert model with single wavelengths for each channel
        SBLR       % regression using Beer-Lambert model with wavelength distributions for each channel
         
                            % structures have fields: 
                            % coef --> 3D array as follows
                            % fluorescence channels are regressed onto
                            % reflectance channels:
                            %   1st dimension sums over each pixel in mask
                            %   2nd dimension sums over each reflectance ch
                            %   3rd dimension sums over each fluorescence ch
                            
                            % labelR --> labels for each reflectance channel
                            % labelF --> labels for each fluorescence channel
                            
                            % varI --> initial variance for each fluor ch pixel before demixing                           
                            % varF --> final variance for each fluor ch pixel after demixing
    end
    
    properties (Constant)
        regTypes = {'BRL','SBLR'}; 
    end
        
    methods
        function obj = BLRClass(cohortName,params)
            obj.cohortName = cohortName;
          
            obj.params = params;
            
            obj.numSessions = 0;
            
            obj.sessionID = cell(0);
            obj.animalID = cell(0);
            obj.folder = cell(0);
            obj.fname = cell(0);    
            
            obj.mask = cell(0);       
            
            obj.BLR = cell(0);
            obj.SBLR = cell(0);
        end
        
        % external Methods
        obj = addSession2BLR(obj,data);
        
    end
end

