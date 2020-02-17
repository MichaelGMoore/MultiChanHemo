function fileNames = scanFolder(folder)
% check that specified folders has .h5 files and the file names

% get the contents of data folder
fdir = dir(folder);

% remove '.' and '..'
fdir = fdir(~ismember({fdir.name},{'.','..'}));

% remove folders
fdir = fdir(~[fdir.isdir]);

% keep only .h5 extensions
fileNames = {};
numFiles = 0;
for f = 1:length(fdir) 
    [filepath name ext] = fileparts(fdir(f).name);
    if isequal(ext,'.h5')
        numFiles = numFiles +1;
        fileNames{end+1,1} = fullfile(folder,fdir(f).name);
    end
end

disp([num2str(numFiles) ' HDF5 imaging session files found'])


end

