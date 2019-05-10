% Script to compute coeffs for Beer-Lambert and Spectral Beer-Lambert
% Models

%% Load the spectral data:

load('opsBLM.mat')

%% 4-wavelength Beer-Lambert Model

method = 1;   % 1 for mean, 2 for maximum likelihood

% mean wavelengths:
if method == 1
    lambda_x = sum(opsBLM.lambdas.*opsBLM.PSex.*opsBLM.PFex)/...
        sum(opsBLM.PSex.*opsBLM.PFex);

    lambda_m = sum(opsBLM.lambdas.*opsBLM.PFem.*opsBLM.PCem)/...
        sum(opsBLM.PFem.*opsBLM.PCem);

    lambda_1 = sum(opsBLM.lambdas.*opsBLM.PS1)/...
        sum(opsBLM.PS1);

    lambda_2 = sum(opsBLM.lambdas.*opsBLM.PS2)/...
        sum(opsBLM.PS2);

% maximum likelihood wavelengths:
elseif method == 2
    [~,idx] = max(opsBLM.PSex.*opsBLM.PFex);
    lambda_x = opsBLM.lambdas(idx);

    [~,idx] = max(opsBLM.PFem.*opsBLM.PCem);
    lambda_m = opsBLM.lambdas(idx);

    [~,idx] = max(opsBLM.PS1);
    lambda_1 = opsBLM.lambdas(idx);

    [~,idx] = max(opsBLM.PS2);
    lambda_2 = opsBLM.lambdas(idx);

end

clear idx method

% the model:

meth = 'makima';
eox = interp1(opsBLM.lambdas,opsBLM.E(:,1),lambda_x,meth,'extrap');
eom = interp1(opsBLM.lambdas,opsBLM.E(:,1),lambda_m,meth,'extrap');
eo1 = interp1(opsBLM.lambdas,opsBLM.E(:,1),lambda_1,meth,'extrap');
eo2 = interp1(opsBLM.lambdas,opsBLM.E(:,1),lambda_2,meth,'extrap');

erx = interp1(opsBLM.lambdas,opsBLM.E(:,2),lambda_x,meth,'extrap');
erm = interp1(opsBLM.lambdas,opsBLM.E(:,2),lambda_m,meth,'extrap');
er1 = interp1(opsBLM.lambdas,opsBLM.E(:,2),lambda_1,meth,'extrap');
er2 = interp1(opsBLM.lambdas,opsBLM.E(:,2),lambda_2,meth,'extrap');

s1m = (eom*er2-erm*eo2)/(eo1*er2-er1*eo2);
s1x = (eox*er2-erx*eo2)/(eo1*er2-er1*eo2);

s2m = (-eom*er1+erm*eo1)/(eo1*er2-er1*eo2);
s2x = (-eox*er1+erx*eo1)/(eo1*er2-er1*eo2);

Xx = .5*interp1(opsBLM.lambdas,opsBLM.x,lambda_x,meth,'extrap');
Xm = .5*interp1(opsBLM.lambdas,opsBLM.x,lambda_m,meth,'extrap');
x1 = interp1(opsBLM.lambdas,opsBLM.x,lambda_1,meth,'extrap');
x2 = interp1(opsBLM.lambdas,opsBLM.x,lambda_2,meth,'extrap');

S_cpwa(1,1) = s1x*(Xx/x1)+s1m*(Xm/x1);
S_cpwa(2,1) = s2x*(Xx/x2)+s2m*(Xm/x2);

clear meth eox eom eo1 eo2 erx erm er1 er2 s1m s1x s2m s2x Xx Xm x1 x2


%% Spectral Beer Lambert Model:

% set the background concentrations:
% Note: the results are not very sensitive to these values, so we just
% choose reasonable estimates
co = (2.2e-3)*.85*.04;
cr = (2.2e-3)*.15*.04;

% set the pathlengths
xm = .5*opsBLM.x;
xx = .5*opsBLM.x;
x1 = opsBLM.x;
x2 = opsBLM.x; 

eo = opsBLM.E(:,1);
er = opsBLM.E(:,2);

ax = sum(opsBLM.PFex.*opsBLM.PSex.*exp(-(co*eo+cr*er).*xx));
daxdco = sum((-eo.*xx).*opsBLM.PFex.*opsBLM.PSex.*exp(-(co*eo+cr*er).*xx));
daxdcr = sum((-er.*xx).*opsBLM.PFex.*opsBLM.PSex.*exp(-(co*eo+cr*er).*xx));

% for emission, we use the Matt Valley 1 camera spectrum
am = sum(opsBLM.PCem.*opsBLM.PFem.*exp(-(co*eo+cr*er).*xm));
damdco = sum((-eo.*xm).*opsBLM.PCem.*opsBLM.PFem.*exp(-(co*eo+cr*er).*xm));
damdcr = sum((-er.*xm).*opsBLM.PCem.*opsBLM.PFem.*exp(-(co*eo+cr*er).*xm));

Ao = damdco/am+daxdco/ax;
Ar = damdcr/am+daxdcr/am;

b1 = sum(opsBLM.PS1.*exp(-(co*eo+cr*er).*x1));
db1dco = sum((-eo.*x1).*opsBLM.PS1.*exp(-(co*eo+cr*er).*x1));
db1dcr = sum((-er.*x1).*opsBLM.PS1.*exp(-(co*eo+cr*er).*x1));
B1o = db1dco/b1;
B1r = db1dcr/b1;


b2 = sum(opsBLM.PS2.*exp(-(co*eo+cr*er).*x2));
db2dco = sum((-eo.*x2).*opsBLM.PS2.*exp(-(co*eo+cr*er).*x2));
db2dcr = sum((-er.*x2).*opsBLM.PS2.*exp(-(co*eo+cr*er).*x2));
B2o = db2dco/b2;
B2r = db2dcr/b2;

denom = B1o*B2r-B1r*B2o;
S_cpa(1,1) = (Ao*B2r-Ar*B2o)/denom;
S_cpa(2,1) = (-Ao*B1r+Ar*B1o)/denom;

clear co cr xm xx x1 x2 eo er ax daxdco daxdcr am damdco damdcr
clear Ao Ar b1 db1dco db1dcr B1o B1r b2 db2dco db2dcr B2o B2r
clear denom

% Results:

% Using Matt Valley 1 Camera emission spectra
% S1 = 1.0626
% S2 = -0.3314

%% Save Results

save('BLM.mat','S_cpa','S_cpwa')





























