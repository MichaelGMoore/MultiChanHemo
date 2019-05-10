function inspectLinearity(obj)
% plot standard deviation and max-value maps for each channel
% Values should be small compared to 1 for linear regression methods to
% apply

fig1 =createFigure(.25,.25,.5,.75);

for k = 1:obj.numCh
    ax1{k} = axes(fig1,'OuterPosition',[(k-1)/obj.numCh,.5,1/obj.numCh,.5]);
    box(ax1{k},'on');
    img = zeros(size(obj.mask));
    img(obj.mask) = max(obj.y{k},[],1);
    imagesc(ax1{k},img,'AlphaData',obj.mask)
    colormap(ax1{k},linspace(0,1,256)'*spectrumRGB(obj.waveCh(k)));   
    colorbar(ax1{k})
    title(ax1{k},{obj.nameCh{k},'Maximum'},'Interpreter','none')
    axis(ax1{k},'image')
    ax2{k} = axes(fig1,'OuterPosition',[(k-1)/obj.numCh,0,1/obj.numCh,.5]);
    box(ax2{k},'on');
    img = zeros(size(obj.mask));
    img(obj.mask) = std(obj.y{k},1,1);
    imagesc(ax2{k},img,'AlphaData',obj.mask)
    colormap(ax2{k},linspace(0,1,256)'*spectrumRGB(obj.waveCh(k)));   
    colorbar(ax2{k})
    title(ax2{k},{obj.nameCh{k},'Standard Deviation'},'Interpreter','none')
    axis(ax2{k},'image')
end
clear k
sgtitle(['Validity of Linearity assumption for ',obj.sessionID,'\_',obj.animalID])

end

