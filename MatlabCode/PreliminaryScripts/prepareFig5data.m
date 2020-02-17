% this code will prepare hdf5 files for the 3 sessions used to generate
% figure5

% we will consider frame-shifts for best alignment

% we provide raw, Beer-Lambert, and MetaModel demixing
0
%% Locate the files
fdir = '/home/mmoore/Optical_Data_Analysis/Data/Widefield/Matt_Valley/GCaMP_Mice/Raw/Sessions4paper';
fname{1} = '171108-M326984_alldat'; % figure 5a
fname{2} = '171018-M287480_alldat'; % figure 5b
fname{3} = '171011-M325951_alldat'; % figure 5c

%% time shifts
startframe{1} = [3,2,1];
startframe{2} = [1,1,1];
startframe{3} = [1,1,1];

%% read each file in and write into new format:
% forgo the spectra

for m = 1:3
    % pull the mouse ID:
    k = strfind(fname{m},'M');
    animalID =  fname{m}(k:(k+6));
    % pull the session ID:
    sessionID = fname{m}(1:6);
    F = h5read(fullfile(fdir,fname{m}),'/JCam/Fluorescence');
    R1 = h5read(fullfile(fdir,fname{m}),'/JCam/Reflectance_575nm');
    R2 = h5read(fullfile(fdir,fname{m}),'/JCam/Reflectance_640nm');  
    % apply time shifts
    F = F(:,:,startframe{m}(1):end);
    R1 = R1(:,:,startframe{m}(2):end);
    R2 = R2(:,:,startframe{m}(3):end);
    % standardize lengths
    N = min([size(F,3),size(R1,3),size(R2,3)]);
    F = F(:,:,1:N);
    R1 = R1(:,:,1:N);
    R2 = R2(:,:,1:N);
    
    % write data to new h5 file
    outFolder = '/home/mmoore/Dropbox/Work/projects/WFOM_MultiChanHemo/Version3/PreliminaryScripts/figure5data';
    outFile = [animalID '.' sessionID '.h5'];
    outName = [outFolder filesep outFile];
    
    h5create(outName,'/images/F',size(F),'ChunkSize',[size(F,1),size(F,2),1],'DataType','uint16');
    h5write(outName,'/images/F',F);
    
    h5create(outName,'/images/R1',size(R1),'ChunkSize',[size(R1,1),size(R1,2),1],'DataType','uint16');
    h5write(outName,'/images/R1',R1);

    h5create(outName,'/images/R2',size(R2),'ChunkSize',[size(R2,1),size(R2,2),1],'DataType','uint16');
    h5write(outName,'/images/R2',R2);
    
    % add some attributes
    
    h5writeatt(outName,'/','animalID',animalID);
    h5writeatt(outName,'/','sessionID',sessionID);
    h5writeatt(outName,'/','CreName','Cux2'); % this will not be a required field
    
end

%% test for time shifts

% f1 = figure(1);
% clf(f1)
% a1 = axes(f1,'outerposition',[0,0,.5,1]);
% a2 = axes(f1,'outerposition',[.5,.5,.5,.5]);
% a3 = axes(f1,'outerposition',[.5 0,.5,.5]);
% 
% img = F(:,:,1);
% imagesc(a1,img)
% axis(a1,'image')
% while true
%     [x,y,button]=ginput(1);
%     if button==3
%         break
%     end
%     x = ceil(x);
%     y = ceil(y);
%     
%     f = squeeze(F(y,x,1:(end-4)));
%     r1 = squeeze(R1(y,x,1:(end-4)));
%     r2 = squeeze(R2(y,x,1:(end-4)));
%     t = (1:length(f))';
%     f = zscore(double(f));
%     r1 = zscore(double(r1));
%     r2 = zscore(double(r2));
%     
%     plot(a2,t,f,'g',t,r1,'y',t,r2,'r')
%     
%     f = squeeze(F(y,x,3:(end-2)));
%     r1 = squeeze(R1(y,x,2:end-3));
%     r2 = squeeze(R2(y,x,1:(end-4)));
%     t = (1:length(f))';
%     f = zscore(double(f));
%     r1 = zscore(double(r1));
%     r2 = zscore(double(r2));
%     
%     plot(a3,t,f,'g',t,r1,'y',t,r2,'r')
%     
% end

