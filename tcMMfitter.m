function tcMMfitter
% fitting of kinetic parameters to experimental results
% includes effects of unstirred layer
% includes ability to read data from Excel spreadsheets
% solute 1: substrate in solution
% solute 2: inhibitor in solution
% solute 3: substrate in cells
% solute 4: inhibitor in cells
% requires Statistics and Machine Learning Toolbox
% TWS, February 2017
% Updated January 2018: near membrane zone removed,
% linear uptake term added, unstirred layer fixed at 0.1 cm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global L DUS DS Jm Kt nconcs concs ntimes times meanconc predconc c0
global funcalls LUS
format compact;
format short e;
% initial estimates of parameters to be fitted
Lcell = 12e-4; %0.01; % Thickness (height) of cell layer in cm
Jl1 = 3.4e-7; % linear (first order) uptake rate in cm/s
Jm1 = 5e-4; % 0.001; % substrate max uptake rate in nmol/cm2/s
Kt1 = 5; % 1; % substrate uptake Michaelis constant (uM = nmol/cm3)

params0 = [Lcell, Jl1, Jm1, Kt1];
% typical values of parameters - needed for setting optimization steps
% use multiple of initial estimates for simplicity
% using a multiple gives more reliable numerical derivatives
refparams = 30.*[Lcell, Jl1, Jm1, Kt1];  
% upper and lower bounds on parameters
ub = [1,1,1,2000];
lb = [0,0,0,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fixed parameters
L = 1; % Total depth of fluid in cm
LUS = 0.1; % Unstirred layer thickness in cm
c01 = 0; % initial substrate concentration in solution (uM = nmol/cm3) (temp. value)
c02 = 0; % initial inhibitor concentration in solution (uM = nmol/cm3)
c03 = 0; % initial substrate concentration in cells (uM = nmol/cm3)
c04 = 0; % initial inhibitor concentration in cells (uM = nmol/cm3)
c0 = [c01 c02 c03 c04];

DUS1 = 6e-6; % substrate diffusivity in cm^2/s
DUS2 = 6e-6; % inhibitor diffusivity in cm^2/s
DS1 = 0.1; % artificial high solute diffusivity in stirred layer, cm^2/s
DS2 = 0.1; % artificial high solute diffusivity in stirred layer, cm^2/s
DS3 = 1.; % artificial high intracellular solute diffusivity, cm^2/s
DS4 = 1.; % artificial high intracellular solute diffusivity, cm^2/s
DUS = [DUS1 DUS2 DS3 DS4];  %no difference in intracellular domain
DS = [DS1 DS2 DS3 DS4];

%Jm1 = 0; % substrate max uptake rate in nmol/cm2/s
%Kt1 = 0; % substrate uptake Michaelis constant (uM = nmol/cm3)
Jm2 = 0; % inhibitor max uptake rate in nmol/cm2/s
Kt2 = 10; % inhibitor uptake Michaelis constant (uM = nmol/cm3)
Jm3 = Jm1; %  0.0012; % substrate max efflux rate in nmol/cm2/s
Kt3 = 5*Kt1; %10; % substrate efflux Michaelis constant (uM = nmol/cm3)
Jm4 = 0.000; % inhibitor max efflux rate in nmol/cm2/s
Kt4 = 250; % inhibitor efflux Michaelis constant (uM = nmol/cm3)
Jm = [Jm1 Jm2 Jm3 Jm4];
Kt = [Kt1 Kt2 Kt3 Kt4];

%read spreadsheet
%if there is a value in the substrate concentration row, data are read!
USLdata = xlsread('Spreadsheet1.xlsx');
size1 = size(USLdata);
ntimes = size1(1) - 2;
ncols = size1(2) - 1;
times = USLdata(3:ntimes + 2,1)';
concsall = USLdata(1,2:ncols + 1);
%nconcs - number of concentrations
%concs - values of concentration
%ndata - number of many data sets for each concentration
nconcs = 0;
concprev = 0;
for icol = 1:ncols
    if concsall(icol) > 0
        if concsall(icol) ~= concprev
            nconcs = nconcs + 1;
            concs(nconcs) = concsall(icol);
            ndata(nconcs) = 1;
            startdata(nconcs) = icol;   %starting column for reading data
            concprev = concsall(icol);
        else
            ndata(nconcs) = ndata(nconcs) + 1;
        end
    end
end
%mean and standard deviation of data, excluding missing values
meanconc = zeros(ntimes,nconcs);
predconc = zeros(ntimes,nconcs);
sdconc = zeros(ntimes,nconcs);
for iconc = 1:nconcs
    for itime = 1:ntimes
        meanconc(itime,iconc) = 0;
        sdconc(itime,iconc) = 0;
        n = 0;
        for icol = 1:ndata(iconc)
            value = USLdata(itime+2,startdata(iconc) + icol);
            if value > 0.
                meanconc(itime,iconc) = meanconc(itime,iconc) + value;
                sdconc(itime,iconc) = sdconc(itime,iconc) + value^2;
                n = n + 1;
            end
        end
        if n > 0
            meanconc(itime,iconc) = meanconc(itime,iconc)/n;
        else
            Disp('*** Error: missing data');
        end
        if n > 1
            sdconc(itime,iconc) = sdconc(itime,iconc)/n;
            sdconc(itime,iconc) = sdconc(itime,iconc) - meanconc(itime,iconc)^2;
            sdconc(itime,iconc) = sdconc(itime,iconc)*n/(n-1);
        else
            sdconc(itime,iconc) = 0;
        end
    end
end

%declare handle for function
mdl = @(params)(Unstirredsolver2017(params));
%run parameter estimation
funcalls = 0;
%The following is for R2016b
 options = optimoptions('lsqnonlin','FiniteDifferenceType','central',...
    'FunctionTolerance',1e-4,'TypicalX',refparams); 
%The following is for R2012a
%options = optimset('FinDiffType','central',...
%    'TolFun',1e-4,'TypicalX',refparams);
paramsest = lsqnonlin(mdl,params0,lb,ub,options)
funcalls
figure(1);
for iconc = 1:nconcs
    plot(times,meanconc(:,iconc),'-o',times,predconc(:,iconc))
    if iconc == 1
        hold on
    end
end
title('Measured and predicted cellular uptake');
xlabel('Time (s)');
ylabel('Uptake (pmol/cm2)');
end

function resid = Unstirredsolver2017(params)
%   Diffusion solver applied to unstirred layer problem
%   2-solute version: substrate and inhibitor
%   Simulates a range of substrate levels with and without inhibitor
%   or a range of inhibitor levels
%  
%   TWS, January 2014. Updated April-May 2014, March 2015, January 2018.
%   This version includes accumulation of substrate in the cells. May 2014.
%   The intracellular concentrations are represented as concentrations in
%   the extracellular domain, but with artificially high diffusion.
%   The cell membrane with Michaelis-Menten uptake is at x = 0
%   In the form expected by PDEPE, the PDE is
%
%       1*  D_ [u] = D_ [D * Du/Dx ] +      0
%           Dt       Dx
%    ---      ---      -------------    -------------
%     c      Du/Dt     f(x,t,u,Du/Dx)   s(x,t,u,Du/Dx)
%
%  The left bc for extracellular is:
%  Net uptake rate + [-1] .* [ D * Du/Dx ] = [0] 
%     ---          ---    -------------     ---
%   p(0,t,u)      q(0,t)   f(0,t,u,Du/Dx)    0
%
%  The left bc for intracellular is:
%  -Net uptake rate + [-1] .* [ D * Du/Dx ] = [0] 
%     ---          ---    -------------     ---
%   p(0,t,u)      q(0,t)   f(0,t,u,Du/Dx)    0
%
%  The right bc for extracellular and intracellular is:
%        0        +  [1] .* [ D * Du/Dx ] = [0] 
%     ---            ---   ---------------   ---
%    p(0,t,u)      q(0,t)    f(0,t,u,Du/Dx)   0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global L Jm Kt nconcs concs ntimes times meanconc predconc LUS c0
global Cellfac funcalls Jl1
%variable parameters
Lcell = params(1); % Thickness of cell layer in cm
Jl1 = params(2); % Linear substrate uptake rate
Jm(1) = params(3); % substrate max uptake rate in nmol/cm2/s
Kt(1) = params(4); % substrate uptake Michaelis constant (uM = nmol/cm3)
Jm(3) = params(3); % substrate max efflux rate in nmol/cm2/s
Kt(3) = 5*params(4); % substrate efflux Michaelis constant (uM = nmol/cm3)
Cellfac = L/Lcell; % Factor to scale for difference in extracellular
% and intracellular compartment sizes
% generate x-vector with close spacing at points of discontinuity
x = [];
if LUS > 0
    xvec = LUS*[0 0.01 0.1 0.3 0.5 0.7 0.9 0.99];
    x = cat(2,x,xvec);
end
if L > LUS
    xvec = LUS + (L - LUS)*[0 0.01 0.1 0.3 0.5 0.7 0.9 0.99 1];
    x = cat(2,x,xvec);
end
%***** start loop over substrate concs ******
t = cat(2,0,times); % has length ntimes + 1
resid(ntimes*nconcs) = 0;   % prellocate memory
rss = 0;
for n = 1:nconcs
    c0(1) = concs(n);
    m = 0;
    sol = pdepe(m,@us1pde,@us1ic,@us1bc,x,t);
    % us = sol(:,:,1); % Extract first component as substrate
    % ui = sol(:,:,2); % Extract second component as inhibitor
    usc = sol(:,:,3); % Extract third component as intracellular substrate
    % uic = sol(:,:,4); % Extract fourth component as intracellular inhibitor 
    for i = 1:ntimes
        % ignore results at t = 0, so use i+1 here
        predconc(i,n) = 1000*trapz(x,usc(i+1,:))/Cellfac;
        rss = rss + (predconc(i,n) - meanconc(i,n))^2;
        resid(i+ntimes*(n-1)) = predconc(i,n) - meanconc(i,n);
    end
end
funcalls = funcalls + 1;
rss
end
% -------------------------------------------------------------------------

function [c,f,s] = us1pde(x,~,~,DuDx)
global DUS LUS DS
  c = [1 1 1 1]';
if x < LUS
  f = DUS'.*DuDx;
else
  f = DS'.*DuDx;
end
  s = [0 0 0 0]';
end

% -------------------------------------------------------------------------

function u0 = us1ic(~)
% initial solute level
global c0
  u0 = c0';
end
    
% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = us1bc(~,ul,~,~,~)
% Michaelis-Menten at x = 0, with inhibition for substrate
% zero flux at x = L
global Jl1 Jm Kt Cellfac
% Jsol = Jm(1)*ul(1)/(Kt(1)*(1 + ul(2)/Kt(2)) + ul(1))... 
%     - Jm(3)*ul(3)/(Kt(3)*(1 + ul(4)/Kt(4)) + ul(3));
% Jinh = Jm(2)*ul(2)/(Kt(2)*(1 + ul(1)/Kt(1)) + ul(2)) ...
%     - Jm(4)*ul(4)/(Kt(4)*(1 + ul(3)/Kt(3)) + ul(4));
%pl = [Jsol Jinh -Jsol*Cellfac -Jinh*Cellfac]';
%version without inhibitor effects
Jsol = Jm(1)*ul(1)/(Kt(1) + ul(1)) - Jm(3)*ul(3)/(Kt(3) + ul(3)) ...
    + Jl1*ul(1);    % include linear uptake term
pl = [Jsol 0 -Jsol*Cellfac 0]';
ql = [-1 -1 -1 -1]';
pr = [0 0 0 0]';
qr = [1 1 1 1]';
end