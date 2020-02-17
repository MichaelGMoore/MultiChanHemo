% read a demixed dataset

%% Read a dataset
fdir = '/home/mmoore/Dropbox/Work/projects/WFOM_MultiChanHemo/Version3/';
fname = 'demixed_noL2.M326984.171108.h5';
inName = [fdir filesep fname];

mask = h5read(inName,'/images/mask');
mask = logical(mask);
F_raw = h5read(inName,'/images/F_raw');
F_MM =  h5read(inName,'/images/F_MM');
F_BL = h5read(inName,'/images/F_BL');

%% Make a comparison figure
x = 87;
y = 83;

time = (1:size(F_raw,1))'/100;

zlist = zeros(size(mask));
zlist(mask) = 1:sum(mask(:));
z = zlist(y,x);

figure(1)
plot(time,F_raw(:,z),'g',time,F_MM(:,z),'k','linewidth',1)
ylabel('\Delta F/F')
xlabel('time [s]')
set(gca,'fontsize',16)
legend('raw','demixed','location','southeast')

clear x y z

%% Deconvolve BOLD
% F_raw = F_MM + T

T_full = F_raw - F_MM;

x = 87;
y = 83;

zlist = zeros(size(mask));
zlist(mask) = 1:sum(mask(:));
z = zlist(y,x);

T = T_full(:,z);
F = F_MM(:,z);

% detrend pixels
T = detrend(T);
F = detrend(F);

% filter
freq = fft_freqs(size(F,1),100);
myfilt =  (abs(freq) < 2) & (abs(freq) > .01);
T = real(ifft(myfilt.*fft(T)));
F = real(ifft(myfilt.*fft(F)));

time = (1:size(F_raw,1))'/100;

figure(1)
plot(time,F,'k',time,-T,'r','linewidth',1)
ylabel('\Delta F/F')
xlabel('time [s]')
set(gca,'fontsize',16)
legend('demixed','Transmittance','location','southeast')

clear x y z

% create a kernel
tbar = 1;
sigma = .5;
k = (tbar/sigma)^2;
theta = sigma^2/tbar;
K = - gampdf(time,k,theta);

figure(2)
plot(time,K,'.-')
xlim([0 10])

% convolve F with K
Tc = conv(K,F,'full');
Tc = Tc(1:length(T));

figure(3)
plot(time,zscore(T),'r',time,zscore(Tc),'b')

%% Fit parameters of gamma function with all pixels
fs = 100;
N = size(F,1);
t = (0:(N-1))'/fs;
k = @(c) -c(1)*gampdf(t,(c(2)/c(3))^2,c(3)^2/c(2));
Tc = @(c) real(ifft(fft(F).*fft(k(c))));

fun = @(c) sum((T - Tc(c)).^2);

c0 = [.1,.5,.25];
cL = [0,0,0];
cU = [1,10,10];
c = fmincon(fun,c0,[],[],[],[],cL,cU);

coef = corrcoef(T,Tc(c));
coef = coef(1,2);

%% Demix using just the gamma function

K = k(c);

I = T + F; % the recorded intensity

% F = I - K*F
% (1+K)*F = I
% F = (1+K)\I;

dom = 1:10000;
col = -K(dom);
col(1) = col(1)+1;
row = zeros(size(col));
row(1) = col(1);
M = toeplitz(col,row);

Fkernel = M\I(dom);

figure
plot(time(dom),I(dom),'g',time(dom),F(dom),'k',time(dom),Fkernel,'m--','linewidth',1)


%% Look for artifacts in original data

addpath(genpath('/home/mmoore/Dropbox/Work/projects/WFOM_MultiChanHemo/Version3'));
fdir = '/home/mmoore/Dropbox/Work/projects/WFOM_MultiChanHemo/Version3/Figures/figure5data';
fname = 'M326984.171108.h5'; % fig5a
name = fullfile(fdir,fname);
fid = H5F.open(name);
gid = H5G.open(fid,'images');
if H5L.exists(gid,'F','H5P_DEFAULT')
    setname = 'F';
else
    setname = 'F1';
end
% read the first fluor channel
I = h5read(fname,['/images/',setname]);
I = permute(I,[3,1,2]);
img = mean(I,1);
img = (img - min(img(:)))/(max(img(:))-min(img(:)));

figure(2)
imagesc(img)
axis image
title('mean image of raw fluorescence')

%% find blood vessels

% manually locate brain
figure(3)
imagesc(img)
axis image
mL = roipoly;
mR = roipoly;

mask = mL | mR;

% find shadows in mean image
blursize = 0.6;
threshold = 0.002;
opensize = 20;
BV = mask.*imgaussfilt(img,blursize)-img > threshold;

% morphological opening to remove the small objects
BV = bwareaopen(BV,opensize,8);

figure(3)
imagesc(img,'alphadata',1-.5*BV)
axis image
title('blood vessel map')

%% compute dF/F

F = double(I);
F = F(:,mask);
F = F./mean(F,1)-1;

%% Inpainting Differential

[X,Y] = meshgrid(1:size(I,3),1:size(I,2));

XnoBV = X(mask & ~BV);
YnoBV = Y(mask & ~BV);
X = X(mask);
Y = Y(mask);

figure(4)
scatter(XnoBV,YnoBV,'b.')
axis image ij
title('Distinguishing Blood Vessel Pixels')
drawnow

FnoBV = F(:,~BV(mask));
Finpaint = zeros(size(F));
% this takes ~1 min:
for f = 1:size(F,1)
    Finpaint(f,:) = griddata(XnoBV,YnoBV,FnoBV(f,:),X,Y);
end
clear f 
ID = F - Finpaint;
ID = ID(:,BV(mask));


%% pull out blood vessel pixels and look at cwt

x = mean(F,2);
fs = 100;
t = (1:size(ID,1))'/fs;

dom = t >= 140& t <= 155;
t = t(dom);
y = x(dom);
freq = fft_freqs(length(x),fs);

wparams = [10,400];
% wparams = [3,60];
fb = cwtfilterbank('SignalLength',length(t),'SamplingFrequency',fs,...
    'FrequencyLimits',[9,15],'VoicesPerOctave',48,'WaveletParameters',wparams);
[s,f,coi] = wt(fb,y);
% [s,f] = wsst(x,fs); % look at wavelet synchro squeezed transform

%%
rowMax = nan(1,size(s,2));
rowf = nan(1,size(s,2));
a = nan(1,size(s,2));
for m = 1:size(s,2)
    [~,i] = max(abs(s(:,m)));
    rowMax(1,m) = i;
    rowf(1,m) = f(rowMax(1,m));
    a(1,m) = s(i,m);
end
clear m



fig = figure(1);
ax=subplot(5,1,[1,2]);
set(ax,'fontsize',16)
XX = [t(1), t(end)];
YY = [min(f),max(f)];
ZZ = zeros(2,2);
surface(ax,XX,YY,ZZ,flipud(abs(s)),...
    'CDataMapping','scaled','FaceColor','texturemap', 'EdgeColor', 'none');
ax.YLim = YY;
ax.XLim = XX;
ax.YDir = 'normal';
ax.YScale = 'log';
xticklabels(ax,{})
yticks(ax,[10,15])
ylabel(ax,'frequency [Hz]')
hold(ax,'on')
plot([t(1),t(end)],13.6*[1 1],'m','linewidth',2)
plot([t(1),t(end)],12.2*[1 1],'m','linewidth',2)
plot([t(1),t(end)],14.0*[1 1],'g','linewidth',2)
plot([t(1),t(end)],10.4*[1 1],'g','linewidth',2)
plot(t,smooth(rowf,15),'c','linewidth',2)
hold(ax,'off')

ax2 = subplot(5,1,3);

filt1 = abs(freq) < 13.6 & abs(freq) > 12.4;
xfilt1 = real(ifft(filt1.*fft(x)));
xfilt1 = xfilt1(dom);
plot(ax2,t,xfilt1,'k','linewidth',.5)
xlim(ax2,[t(1),t(end)])
ylim(ax2,[-.003 .003])
xticklabels(ax2,{})
yticklabels(ax2,{})
set(ax2,'fontsize',16)
box(ax2,'on')

ax3 = subplot(5,1,4);
filt2 = abs(freq) < 14.0 & abs(freq) > 10.4;
xfilt2 = real(ifft(filt2.*fft(x)));
xfilt2 = xfilt2(dom);
plot(ax3,t,xfilt2,'k','linewidth',.5)
xlim(ax3,[t(1),t(end)])
ylim(ax3,[-.003 .003])
xticklabels(ax3,{})
yticklabels(ax3,{})
set(ax3,'fontsize',16)
box(ax3,'on')

ax4 = subplot(5,1,5);
plot(t-t(1),real(a),'k','linewidth',.5)
xlim(ax4,[0,t(end)-t(1)])
ylim(ax4,[-.003 .003])
xlabel(ax4,'time [s]')
yticklabels(ax4,{})
set(ax4,'fontsize',16)
box(ax4,'on')

