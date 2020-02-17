% Prepare Rorb GFP data as .h5 with simple schema
% Schema: Reflectance Data + Fluorescence data + spectra

addpath(genpath('/home/mmoore/Dropbox/Work/projects/WFOM_MultiChanHemo/Version3'));

folder = '/home/mmoore/Optical_Data_Analysis/Data/Widefield/Matt_Valley/GFP_Mice/Rorb/Raw';
fdir = dir(folder);
fdir = fdir(3:end);

% read each .h5 file and write a new .h5 using simplified schema
load('mySpectra.mat')
for m = 1:length(fdir)
    info = h5info([folder filesep fdir(m).name]);
    % pull the mouse ID:
    k = strfind(fdir(m).name,'M');
    animalID =  fdir(m).name(k:(k+6));
    % pull the session ID:
    sessionID = fdir(m).name(1:6);
    % pull the fluoresence data
    F = h5read([folder,filesep,fdir(m).name],'/JCam/Fluorescence');
    R1 = h5read([folder,filesep,fdir(m).name],'/JCam/Reflectance_575nm');
    R2 = h5read([folder,filesep,fdir(m).name],'/JCam/Reflectance_640nm');

    % write data to new h5 file

    outFolder = '/home/mmoore/Dropbox/Work/projects/WFOM_MultiChanHemo/Version3/PreliminaryScripts/Rorb';
    outFile = [animalID '.' sessionID '.h5'];
    fname = [outFolder filesep outFile];
    
    h5create(fname,'/images/F',size(F),'ChunkSize',[size(F,1),size(F,2),1],'DataType','uint16');
    h5write(fname,'/images/F',F);
    
    h5create(fname,'/images/R1',size(R1),'ChunkSize',[size(R1,1),size(R1,2),1],'DataType','uint16');
    h5write(fname,'/images/R1',R1);

    h5create(fname,'/images/R2',size(R2),'ChunkSize',[size(R2,1),size(R2,2),1],'DataType','uint16');
    h5write(fname,'/images/R2',R2);
    
    % add some attributes
    
    h5writeatt(fname,'/','animalID',animalID);
    h5writeatt(fname,'/','sessionID',sessionID);
    h5writeatt(fname,'/','CreName','Rorb'); % this is not be a required field
    
    % Add spectra data
    
    h5create(fname,'/spectra/Fex',size(spectra.Fex),'DataType','double');
    h5write(fname,'/spectra/Fex',spectra.Fex);
    h5create(fname,'/spectra/Fem',size(spectra.Fem),'DataType','double');
    h5write(fname,'/spectra/Fem',spectra.Fem);    
    h5create(fname,'/spectra/R1',size(spectra.R1),'DataType','double');
    h5write(fname,'/spectra/R1',spectra.R1);      
    h5create(fname,'/spectra/R2',size(spectra.R2),'DataType','double');
    h5write(fname,'/spectra/R2',spectra.R2);
    
end
clear m k
clear fdir fname folder outFolder outFile
clear animalID sessionID info spectra F R1 R2
    
    


