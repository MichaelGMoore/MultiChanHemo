classdef MMRClass
    % This class contains Meta-Model regression data for a cohert of GFP
    % reporter mice
    
    properties
        cohortName          % a name for the PWR cohort, e.g. 'Cux2'       

        params              % all pipeline parameters used in processing the data
        
        numSessions         % number of sessions included in the PWR cohort
                            %   updated dynamically as each session is added                                  
              
        trainedMM           % meta-model coefficients trained on all data submitted
        
        % cell array properties (cell labels each session)
        sessionID           % ID for each session
        animalID            % ID for each animal
        folder              % folder containing raw data
        fname               % raw data filename
        
        mask                % mask for each animal 
                            %   D = number of pixels in mask
        
        % meta-model training data
        MMdata              % meta-model data for each session
                                                        
        % demixing results
        MMR                 % sessions demixed on their own metamodels
        CVMMR               % sessions demixed using metamodels trained on other animals only
               
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
        regTypes = {'MMR','CVMMR'};
        
    end
       
    methods
        function obj = MMRClass(cohortName,params)
            obj.cohortName = cohortName;
          
            obj.params = params;
            
            obj.numSessions = 0;
            
            obj.trainedMM = [];
            
            obj.sessionID = cell(0);
            obj.animalID = cell(0);
            obj.folder = cell(0);
            obj.fname = cell(0);    
            
            obj.mask = cell(0);     
           
            obj.MMdata = cell(0);
            
            obj.MMR = cell(0);
            obj.CVMMR = cell(0);
      
        end
          
        % external Methods
        obj = addSession2MMR(obj,data);
        
        % once all sessions are added, these can be run
        obj = computeMMR(obj,data);
        obj = computeCVMMR(obj,data);       
        obj = trainMM(obj);
    end
end

