function inspectDemixing(obj,data,type)
assert(isa(data,'dataClass'),'data must be dataClass object')


% find the session
for m = 1:length(obj.sessionID)
    sesID = obj.sessionID{m};
    anID = obj.animalID{m};
    if isequal(data.sessionID,sesID) && isequal(data.animalID,anID)
        break
    end
end
% now m is the session index

% compute demixed signal
fluoCh = find(contains(data.typeCh,'F'));
reflCh = find(contains(data.typeCh,'R'));
demixed = zeros(data.frames,sum(data.mask(:)),length(fluoCh));
S = getfield(obj,type);
% S = S{m};

for f = 1:length(fluoCh)
    raw = data.y{fluoCh(f)};
    approx = zeros(size(raw));
    for r = 1:length(reflCh)
        approx = approx + data.y{reflCh(r)}.*S.coef(:,r,f)';
    end
    demixed(:,:,f) = raw - approx;
end

% compare raw and demixed
N = data.frames;
D = sum(data.mask(:));

f1 = createFigure(.1,.1,.8,.8);
ax1 = axes(f1,'outerposition',[0,0,.5,1],'box','on');
for f = 1:length(fluoCh)
    ax2{f} = axes(f1,'outerposition',[.5,(f-1)/length(fluoCh),.5,1/length(fluoCh)],'box','on');
end

img = zeros(size(data.mask));
img(data.mask) = S.varI(:,1);
imagesc(ax1,img,'AlphaData',data.mask)
axis(ax1,'image')
set(ax1,'CLim',[0,prctile(img(:),99.9)])
title(ax1,{[data.sessionID,'\_',data.animalID],['Initial variance of ',S.labelF{1}]})
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
            Y = raw(:,d);
            plot(ax2{f},data.times,Y,'Color',spectrumRGB(data.waveCh(fluoCh(f))),'LineWidth',1)
            Y = demixed(:,d);
            plot(ax2{f},data.times,Y,'k','LineWidth',1)
            hold(ax2{f},'off')
        end

        drawnow
    end
    
end


end

