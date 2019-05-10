function [images,mask,metaData] = readMattValleyh5(folder,fname,frames2read)
% Custom Parser for MattValley h5 3 channel datafiles
% Inputs:
%   folder          folder name with no trailing '/'
%   name            file name (including any extension)
%   frames2read     array of frame indices to include 
%                   (with 2 input variables will read all frames)

% Outputs:
%   images      data cell array. data{ch} will be data for channel in (t,y,x) form
%   mask        mask of brain region
%   metaData    structure containing metadata 
%       fields:
%           folder      folder where data is located
%           fname       filename of raw data
%           anID        animal id
%           sessionID   session id
%           numCh       number of channels
%           nameCh      [numCh x 1] cell array name of channels
%           waveCh      [numCh x 1] array of approximate channel wavelengths
%           typeCh      [numCh x 1] cell array, 'F' or 'R' to indicate Fluorescence or Reflectance
%           fs          sampling frequency
metaData = struct;
metaData.folder = folder;
metaData.fname = fname;

k = strfind(fname,'M');
metaData.animalID = fname(k:(k+6)); % M-number mouse label
metaData.sessionID = fname(1:6); % date of session

info = h5info([folder,filesep,fname]);

metaData.numCh = length(info.Groups(1).Datasets);
for k = 1:metaData.numCh
    metaData.nameCh{k} = info.Groups(1).Datasets(k).Name;    
    images{k} = h5read([folder,filesep,fname],[info.Groups(1).Name,'/',info.Groups(1).Datasets(k).Name]);
    if nargin == 3
        images{k} = images{k}(:,:,frames2read);
    end
    images{k} = permute(double(images{k}),[3,1,2]);
end
metaData.waveCh = [520,575,640];
metaData.typeCh = {'F','R','R'};
metaData.fs = 100;

% initial mask based on mean intensity of fluorescence channel
img = mean(images{1},1);
img = squeeze(img);
img = (img - min(img(:)))/(max(img(:))-min(img(:)));
mask = img > .1;
mask = bwareaopen(mask,2000);
mask = imfill(mask,'holes');
mask = bwmorph(mask,'open');

% update mask to remove outlier pixels
% for ch = 1:metaData.numCh
%     % for each channel compute std
%     sigma{ch} = squeeze(std(images{ch},1,1)./mean(std(images{ch},1)));
%     % for each channel compute the median std
%     medianSigma{ch} = median(sigma{ch}(:));
%     % flag all pixels about 10*medianSigma
%     flagCh{ch} = sigma{ch} > 5*medianSigma{ch};
% end
% % merge the channels
% flag = false(size(sigma{1}));
% for ch = 1:metaData.numCh
%     flag = flag | flagCh{ch};
% end
% % dilate the flagged regions by 1 pixel
% % flag = bwmorph(flag,'thicken',1);
% % update the mask
% mask = mask & ~flag;

% createFigure(.25,.1,.5,.8);
% imagesc(img,'AlphaData',.8+.2*mask)
% axis image
% title(fname,'Interpreter','none')
% drawnow


end













