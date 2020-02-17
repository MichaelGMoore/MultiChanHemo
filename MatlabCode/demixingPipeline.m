% Pipeline to demix a dataset using pre-trained models

% add path to project
addpath(genpath('/home/mmoore/Dropbox/Work/projects/WFOM_MultiChanHemo/Version3'));

% specify the file to be demixed
fdir = '/home/mmoore/Dropbox/Work/projects/WFOM_MultiChanHemo/Version3/Figures/figure5data';
fname{1} = 'M326984.171108.h5'; % fig5a
fname{2} = 'M287480.171018.h5'; % fig5b
fname{3} = 'M325951.171011.h5'; % fig5c

%% use the Cux2 trained models for demixing
load('SCux2'); % load the datastructure for the training data
S = SCux2;

for f = 1
    % Read the data
    name = fullfile(fdir,fname{f});
    [y,mask,stats,meanFImg,meanRImg] = processH5data(name,S);

    figure(f)
    colormap plasma
    img = zeros(size(mask));
    img(mask) = std(y(:,:,1),1,1);
    imagesc(img,'alphadata',mask)
    axis image

    % compute metamodel predictors
    predictors = computeMMP(stats,meanFImg,meanRImg,mask,S);

    % Demix with Beer-Lambert
    C_BL = S.session(1).BL.coef(1,:);
    F_BL = y(:,:,1) - C_BL(1)*y(:,:,2) - C_BL(2)*y(:,:,3);
    
    % rename raw data
    F_raw = y(:,:,1);

    % Demix with Meta-Model
    P = zscore(predictors,1,1);

    C = S.CVMM0.ResponseMean + P*S.CVMM0.MMcoef;

    F_MM = zeros(size(y(:,:,1)));
    for d = 1:sum(mask(:))
        F_MM(:,d) = y(:,d,1)-C(d,1)*y(:,d,2)-C(d,2)*y(:,d,3);
    end
    
    % write to file
    animalID = fname{f}(1:7);
    sessionID = fname{f}(9:14);
    
    outFolder = '/home/mmoore/Dropbox/Work/projects/WFOM_MultiChanHemo/Version3/';
    outFile = ['demixed_noL2.',animalID '.' sessionID '.h5'];
    outName = [outFolder filesep outFile];
    
    h5create(outName,'/images/mask',size(mask),'DataType','double');
    h5write(outName,'/images/mask',double(mask));   
    
    h5create(outName,'/images/F_raw',size(F_raw),'ChunkSize',[1,size(F_raw,2)],'DataType','double');
    h5write(outName,'/images/F_raw',F_raw);
    
    h5create(outName,'/images/F_BL',size(F_BL),'ChunkSize',[1,size(F_BL,2)],'DataType','double');
    h5write(outName,'/images/F_BL',F_BL);
    
    h5create(outName,'/images/F_MM',size(F_MM),'ChunkSize',[1,size(F_MM,2)],'DataType','double');
    h5write(outName,'/images/F_MM',F_MM);

end

%% Examine the pixels

f2 = figure(2);
colormap plasma
a1 = axes(f2,'outerposition',[0,0,.5,1]);
a2 = axes(f2,'outerposition',[.5,.5,.5,.5]);
a3 = axes(f2,'outerposition',[.5,0,.5,.5]);

title(a2,'Beel-Lambert demixing')
set(a2,'fontsize',16)
title(a3,'Meta-Model demixing') 
set(a3,'fontsize',16)

img = zeros(size(mask));
img(mask) = std(F_raw,1,1);
imagesc(a1,img,'alphadata',mask)
axis(a1,'image')

zlist = zeros(size(mask));
zlist(mask) = 1:sum(mask(:));
while true
    [xx,yy,bb] = ginput(1);
    if bb == 3
        break
    end
    xx = ceil(xx);
    yy = ceil(yy);
    if mask(yy,xx)
        plot(a2,F_raw(:,zlist(yy,xx)),'g')
        hold(a2,'on')
        plot(a2,F_BL(:,zlist(yy,xx)),'k')
        hold(a2,'off')
        title(a2,'Beel-Lambert demixing')
        
        plot(a3,F_raw(:,zlist(yy,xx)),'g')
        hold(a3,'on')
        plot(a3,F_MM(:,zlist(yy,xx)),'k')
        hold(a3,'off')
        title(a3,'Meta-Model demixing') 
              
    end
end






