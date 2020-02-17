function [y,mask,stats,meanFImg,meanRImg] = processH5data(fname,S)
% 10/30/2019
% Michael G. Moore, Michigan State University

% take an HDF5 file in the proper format, and do the pre-processing
% Preprossing steps:
% (1) load datameanFImg
% (2) make a mask
% (3) compute and save mean images
% (4) compute y = dI/I for each channel
% (5) apply low-pass filter if specified by project

% check if first fluor channel is 'F' or 'F1'
fid = H5F.open(fname);
gid = H5G.open(fid,'images');
if H5L.exists(gid,'F','H5P_DEFAULT')
    setname = 'F';
else
    setname = 'F1';
end
% read the first fluor channel
I = h5read(fname,['/images/',setname]);
I = permute(I,[3,1,2]);
% construct the mean image for blood vessel mapping
img = squeeze(mean(I,1));

% make a mask of brain-region
mask = makeMask(img);
meanFImg{1} = img;

% load all channels into a y(time,pixel,channel) array
%   pre-allocate
y = zeros(size(I,1),sum(mask(:)),S.numF+S.numR);
% insert the first channel
y(:,:,1) = I(:,mask);
% read the remaining fluor channels
for f = 2:S.numF
    setname = ['F',num2str(f)];
    I = h5read(fname,['/images/',setname]);
    I = permute(I,[3,1,2]);
    img = squeeze(mean(I,1));
    meanFImg{f} = img;
    y(:,:,f) = I(:,mask);
end
% read the reflectance channels
for r = 1:S.numR
    setname = ['R',num2str(r)];
    I = h5read(fname,['/images/',setname]);
    I = permute(I,[3,1,2]);
    img = squeeze(mean(I,1));
    meanRImg{r} = img;
    y(:,:,S.numF+r) = I(:,mask);
end

% compute dI/I
y = y./mean(y,1)-1;

% apply low pass filter
if S.params.LPF.doLPF
    n = 3;
    nyq = S.Fs/2;
    Wn = S.params.LPF.f0/nyq;
    [b,a] = butter(3,Wn);
    y = filtfilt(b,a,y);
end

% %% find 10-sigma outliers in varI and remove pixels within 2 pixels of the outlier
% 
% varI = squeeze(var(y,0,1));
% bad = any(isoutlier(varI,'mean',1,'ThresholdFactor',10),2); 
% % make an image of the bad pixels
% badMask = zeros(size(mask),'logical');
% badMask(mask) = bad;
% % dilate image
% badMask = imdilate(badMask,strel('disk',2));
% % convert back to bad array
% bad = badMask(mask);  
% 
% y = y(:,~bad,:);
% mask = mask & ~badMask;

% collect pixel-wise statistics for meta-model
N = size(y,1);
D = size(y,2);

R = y(:,:,1:S.numF);           % [N x D x numF], Responses = fluo channels
P = y(:,:,S.numF+(1:S.numR));  % [N x D x numR], Predictors = reflectance channels

stats = struct;
stats.L1 = squeeze(mean(abs(P),1)); % [D x numR], L1 norm of Predictors
stats.P2 = squeeze(mean(P.^2,1)); % [D x numR], mean square of Predictors
stats.P3 = squeeze(mean(P.^3,1)); % [D x numR], mean cube of Predictors
stats.P4 = squeeze(mean(P.^4,1)); % [D x numR], mean quart of Predictors

stats.PTP = permute(sum(permute(P,[3,1,2,4]).*permute(P,[4,1,2,3]),2),[1,3,4,2])/N; 
    % [numR x D x numR] pixel-wise correlation matrix of predictors
stats.PTR = permute(sum(permute(P,[3,1,2,4]).*permute(R,[4,1,2,3]),2),[1,3,4,2])/N; 
    % [numR x D x numF] pixel-wise cross correlation matrix of
    % predictors and responses
stats.RTR = permute(sum(permute(R,[3,1,2,4]).*permute(R,[4,1,2,3]),2),[1,3,4,2])/N; 
    % [numF x D x numF] pixel-wise correlation matrix of responses

end

