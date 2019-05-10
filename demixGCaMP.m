% script to use trained Meta Model to demix a GCaMP dataset

% folder = '/home/mmoore/Optical_Data_Analysis/Data/Widefield/Matt_Valley/GCaMP_Mice/Raw/Sessions4paper';
% fname = 'M377518_128x128_alldat';

folder = '/home/mmoore/Optical_Data_Analysis/Data/Widefield/Matt_Valley/GCaMP_Mice/Raw/Sessions4paper';
fname = '171108-M326984_alldat';

[images,mask,metaData] = readMattValleyh5(folder,fname,1:59937);

%% convert to dataClass object

clear data
data = dataClass(images,mask,metaData);
clear images mask metaData

%% Copy the data before filtering

dataRaw = data;

%% Low Pass filter the data

params = struct;
% Low Pass Filter Parameters
params.filter.f0 = 3;  % filter frequency in Hz
params.filter.width = 0; % width of the filter edge in Hz
data = applyLowPassFilter(data,params);

%% Find the best timing alignment
% positive numbers move the peaks to the right
timeShifts = [0,0,0];
inspectTimeShifts(data,timeShifts)

%% Preprocess the data
% if you find timeShifts that you believe are more plausible, enter them below to apply them
params.timeShifts = [0,0,0];
for m = 1:length(data.y)
    data.y{m} = circshift(data.y{m},-params.timeShifts(m),1);
    dataRaw.y{m} = circshift(dataRaw.y{m},params.timeShifts(m),1);
end
clear m

%% Final look at alignment
inspectData(dataRaw)

%% construct the metamodel
% we don't need a new  metamodel, but we need the metamodel predictors

% MMR params
params.MMR.numScales = 6; % number of scales for blood vessel maps
params.MMR.FVEThresh = .8; % threshold on PWR FVE for inclusing in MM training data

MMRGCaMP = MMRClass('Cux2',params); % Meta-Model Regression
MMRGCaMP= addSession2MMR(MMRGCaMP,data);

%% Load the pre-trained meta model
load(['SavedObjects',filesep,'MMRCux2.mat'])

%% demix the data

D = sum(data.mask(:));
N = data.frames;

FVEThresh = params.MMR.FVEThresh;

% parse the reflectance and fluorescence
reflCh = find(contains(data.typeCh,'R'));
fluoCh = find(contains(data.typeCh,'F'));
numReflCh = length(reflCh);
numFluoCh = length(fluoCh);
labelR = cell(numReflCh,1); 
labelF = cell(1,numFluoCh); 
for nr = 1:numReflCh
    labelR{nr,1} = [num2str(data.waveCh(reflCh(nr))),'nm'];
end
clear nr
for nf = 1:numFluoCh
    labelF{1,nf} = [num2str(data.waveCh(fluoCh(nf))),'nm'];
end
clear nf

F = zeros(N,D,numFluoCh);
varI = zeros(D,numFluoCh);
varF = zeros(D,numFluoCh);
FVE = zeros(D,numFluoCh);
for nf = 1:numFluoCh 
    RespMean = MMR.trainedMM.RespMean;
    Coeffs = MMR.trainedMM.Coeffs;
%     RespMean = MMR.CVMMR{1}.RespMean;
%     Coeffs = MMR.CVMMR{1}.Coeffs;
    
    Predictors = MMRGCaMP.MMdata{1}.predictors; % for demixing we need to use all pixels
    
    coef = RespMean + Coeffs*zscore(Predictors,1,2);
    clear RespMean Coeffs Predictors
    
    f = fluoCh(nf);
    for d = 1:D
        R = [];
        for r = reflCh
            R = cat(2,R,data.y{r}(:,d));
        end
        clear r
        FF = dataRaw.y{f}(:,d);
        varI(d,nf) = var(FF,1,1);
        varF(d,nf) = var(FF-R*coef(:,d),1,1);
        FVE(d,nf) = 1-varF(d,nf)/varI(d,nf); 
        F(:,d,nf) = FF-R*coef(:,d);
    end
    clear d R FF 
end
clear nf f

%% Spectral Beer Lambert Regression
% for comparison

params.BLR.SBLRcoefs = [1.176023650535313;-0.434043363539705];

F_SBLR = zeros(N,D,numFluoCh);
for nf = 1:numFluoCh 

    coef = params.BLR.SBLRcoefs;
    
    f = fluoCh(nf);
    for d = 1:D
        R = [];
        for r = reflCh
            R = cat(2,R,data.y{r}(:,d));
        end
        clear r
        FF = dataRaw.y{f}(:,d);
        F_SBLR(:,d,nf) = FF-R*coef(:);
    end
    clear d FF R
end
clear nf f coef

%% Look at demixing:

f1 = createFigure(.1,.1,.8,.8);
ax1 = axes(f1,'outerposition',[0,0,.5,1],'box','on');
for f = 1:length(fluoCh)
    ax2{f} = axes(f1,'outerposition',[.5,(f-1)/length(fluoCh),.5,1/length(fluoCh)],'box','on');
end

img = zeros(size(data.mask));
img(data.mask) = varI(:,1);
imagesc(ax1,img,'AlphaData',data.mask)
axis(ax1,'image')
set(ax1,'CLim',[0,prctile(img(:),99.9)])
title(ax1,{[data.sessionID,'\_',data.animalID],['Initial variance of ',labelF{1}]})
colormap(ax1,linspace(0,1,256)'*spectrumRGB(data.waveCh(fluoCh(1))))
colorbar(ax1)
    
dmap = zeros(size(data.mask));
dmap(data.mask) = 1:D;
while true
    [x,y,button] = ginput(1);
    if button==3
        break
    end
    x = round(x);
    y = round(y);
        
    if data.mask(y,x)

        d = dmap(y,x);

        for f = 1:length(fluoCh)
            cla(ax2{f})
            hold(ax2{f},'on')
            Y = data.y{1}(:,d);
            plot(ax2{f},data.times,Y,'Color',spectrumRGB(data.waveCh(fluoCh(f))),'LineWidth',1)
            Y = F_SBLR(:,d);
            plot(ax2{f},data.times,Y,'r','LineWidth',1)
            Y = F(:,d);
            plot(ax2{f},data.times,Y,'k','LineWidth',1)
            hold(ax2{f},'off')
        end

        drawnow
    end
    
end


%% H5 write

mask = double(data.mask);
sessionID = data.sessionID;
animalID = data.animalID;
filtered = data.filtered;
timeShifts = params.timeShifts;
G = data.y{1};

outfile = ['H5Objects',filesep,sessionID,'_',animalID,'_MM_GCaMP_LargeMidlineArtifact'];

h5create(outfile,'/Data/filtered',size(filtered),'ChunkSize',[1,size(filtered,2)],'DataType','double');
h5write(outfile,'/Data/filtered',filtered)

h5create(outfile,'/Data/timeShifts',size(timeShifts),'ChunkSize',[1,size(timeShifts,2)],'DataType','double');
h5write(outfile,'/Data/timeShifts',timeShifts)

h5create(outfile,'/Data/mask',size(mask),'ChunkSize',[1,size(mask,2)],'DataType','double');
h5write(outfile,'/Data/mask',mask)

h5create(outfile,'/Data/G',size(G),'ChunkSize',[1,size(G,2)],'DataType','double');
h5write(outfile,'/Data/G',G)

h5create(outfile,'/Data/F',size(F),'ChunkSize',[1,size(F,2)],'DataType','double');
h5write(outfile,'/Data/F',F)

h5create(outfile,'/Data/F_SBLR',size(F_SBLR),'ChunkSize',[1,size(F_SBLR,2)],'DataType','double');
h5write(outfile,'/Data/F_SBLR',F_SBLR)
