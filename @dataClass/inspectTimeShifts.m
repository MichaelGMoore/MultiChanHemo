function inspectTimeShifts(obj,timeShifts)

N = obj.frames;
D = sum(obj.mask(:));

% interactive figure to look at data
f1 = createFigure(.1,.1,.8,.8);

ax1 = axes(f1,'outerposition',[0,0,.5,1],'box','on');
ax2 = axes(f1,'outerposition',[.5,0,.5,.5],'box','on');
ax3 = axes(f1,'outerposition',[.5,.5,.5,.5],'box','on');

% find the fluorescence channel
cF = getFluoChan(obj);
cF = cF(1);
colorF = spectrumRGB(obj.waveCh(cF));
cmF = linspace(0,1,256)'*colorF;


img = zeros(size(obj.mask));
img(obj.mask) = std(obj.y{cF},1,1);
imagesc(ax1,img,'AlphaData',obj.mask)
axis(ax1,'image')
set(ax1,'CLim',[0,prctile(img(:),99.9)])
title(ax1,'STD of Fluorescence Channel')
colormap(ax1,cmF)
colorbar(ax1)

hold(ax2,'on')
for c = 1:obj.numCh
    Y = mean(obj.y{c},2);
    plot(ax2,obj.times,Y,'.-','Color',spectrumRGB(obj.waveCh(c)),'LineWidth',1)
end
clear c
hold(ax2,'off')

freqs = fft_freqs(size(obj.y{1},1),obj.fs);

hold(ax3,'on')
for c = 1:obj.numCh
    Y = mean(obj.y{c},2);
    plot(ax3,freqs,abs(fft(Y)),'.-','Color',spectrumRGB(obj.waveCh(c)),'LineWidth',1)
end
clear c
xlim(ax3,[.5,obj.fs/2])
hold(ax3,'off')
title(ax3,'averaged over pixels')
drawnow

dmap = zeros(size(obj.mask));
dmap(obj.mask) = 1:D;

while true
    [x,y,button] = ginput(1);
    if button==3
        break
    end
    x = round(x);
    y = round(y);
        
    if obj.mask(y,x)

        d = dmap(y,x);
        cla(ax2)
        cla(ax3)
        hold(ax2,'on')
        for c = 1:obj.numCh
            Y = circshift(obj.y{c}(:,d),timeShifts(c),1);
            plot(ax2,obj.times,Y,'.-','Color',spectrumRGB(obj.waveCh(c)),'LineWidth',1)
        end
        clear c
        hold(ax2,'off')
        

        freqs = fft_freqs(size(obj.y{1},1),obj.fs);

        hold(ax3,'on')
        for c = 1:obj.numCh
            Y = obj.y{c}(:,d);
            plot(ax3,freqs,abs(fft(Y)),'.-','Color',spectrumRGB(obj.waveCh(c)),'LineWidth',1)
        end
        clear c
        xlim(ax3,[.5,obj.fs/2])
        hold(ax3,'off')   
        
        title(ax3,['y = ',num2str(y),', x = ',num2str(x)])
        drawnow
    end
    
end


end

