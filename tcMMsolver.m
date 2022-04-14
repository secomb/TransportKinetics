function tcMMsolver
%   Diffusion solver applied to unstirred layer problem
%   2-solute version: substrate and inhibitor
%   Simulates net substrate/inhibitor accumulation over time
%   for a range of substrate levels with and without inhibitor
%   or a range of inhibitor levels
%   Includes solute-free zone in initial condition
%   and the impact of mediated efflux of accumulated substrate on net substrate accumulation
%   Also allows ramped or smooth solute concentration within the unstirred layer
%   TWS, January 2014. Updated April-May 2014, March 2015.
%   Updated January 2018 to include a linear (first order) uptake term.
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

global L LUS c0 DUS DS Jl1 Jm Kt L00 ramp Cellfac
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varyinhib = 1;
% varyinhib = 0 runs a range of substrate concentrations,
% with a fixed inhibitor concentration
% varyinhib = 1 runs a fixed substrate concentration,
% with a range of inhibitor concentrations
conc_fac = [1];
%conc_fac = [0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100 300];
% Runs are done for each of the above factors multiplying the substrate
% or inhibitor initial concentration (c01 or c02)
% *** Do not set conc_fac to zero ***
nconcs = length(conc_fac); % Length of conc_fac defines number of runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 1.0; % Total depth of fluid over cell layer in cm
Lcell = 0.0012; % Thickness (height) of cell layer in cm
Cellfac = L/Lcell; % Factor to scale for difference in extracellular
% and intracellular compartment sizes
L00 = 0.0;  % Thickness of initial solute-free/depleted NMZ (near membrane zone)
ramp = 0; % ramp = 0 for zero conc., ramp = 1 for ramped conc. in  NMZ
% ramp = 2, 3 etc. for smooth (power-law) initial profile in NMZ
LUS = 0.100; % Unstirred layer thickness in cm
c01 = 1.0; % initial substrate concentration in solution (uM = nmol/cm3)
c02 = 0; % initial inhibitor concentration in solution (uM = nmol/cm3)
c03 = 0; % initial substrate concentration in cells (uM = nmol/cm3)
c04 = 0; % initial inhibitor concentration in cells (uM = nmol/cm3)

DUS1 = 6e-6; % substrate diffusivity in cm^2/s
DUS2 = 6e-6; % inhibitor diffusivity in cm^2/s
DS1 = 0.1; % artificial high solute diffusivity in stirred layer, cm^2/s
DS2 = 0.1; % artificial high solute diffusivity in stirred layer, cm^2/s
DS3 = 1.; % artificial high intracellular solute diffusivity, cm^2/s
DS4 = 1.; % artificial high intracellular solute diffusivity, cm^2/s

Jl1 = 0; % linear uptake rate in cm/s - added January 2018
Jm1 = 5e-4; % substrate max uptake rate in nmol/cm2/s
Kt1 = 5; % substrate uptake Michaelis constant (uM = nmol/cm3)
Jm2 = 0; % inhibitor max uptake rate in nmol/cm2/s
Kt2 = 10; % inhibitor uptake Michaelis constant (uM = nmol/cm3)
Jm3 = Jm1; % substrate max efflux rate in nmol/cm2/s
Kt3 = 5*Kt1; % substrate efflux Michaelis constant (uM = nmol/cm3)
Jm4 = 0.0; % inhibitor max efflux rate in nmol/cm2/s
Kt4 = 50; % inhibitor efflux Michaelis constant (uM = nmol/cm3)

Tmax = 180 ; % simulation time in s

DUS = [DUS1 DUS2 DS3 DS4];  %no difference in intracellular domain
DS = [DS1 DS2 DS3 DS4];
Jm = [Jm1 Jm2 Jm3 Jm4];
Kt = [Kt1 Kt2 Kt3 Kt4];
c0 = [c01 c02 c03 c04];

% generate x-vector with close spacing at points of discontinuity
Lint1 = min(L00,LUS);
Lint2 = max(L00,LUS);
x = [];
if Lint1 > 0
    xvec = Lint1*[0 0.01 0.1 0.3 0.5 0.7 0.9 0.99];
    x = cat(2,x,xvec);
end
if Lint2 > Lint1
    xvec = Lint1 + (Lint2 - Lint1)...
        *[0 0.01 0.03 0.05 0.07 0.1 0.15 0.2 0.25 0.3 0.5 0.7 0.9 0.99];
    x = cat(2,x,xvec);
end
if L > Lint2
    xvec = Lint2 + (L - Lint2)*[0 0.01 0.1 0.3 0.5 0.7 0.9 0.99 1];
    x = cat(2,x,xvec);
end

t = Tmax*[0 0.001 0.002 0.005 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1...
    0.12 0.14 0.16 0.18 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
ttt = length(t);

format compact;

% prellocate memory
sub_conc(nconcs) = 0; 
inh_conc(nconcs) = 0;
sub_uptake_rate(nconcs) = 0;
inh_uptake_rate(nconcs) = 0;
cell_sub_content(ttt) = 0;
cell_inh_content(ttt) = 0;

%***** start loop over substrate or inhibitor concentrations ******
for n = 1:nconcs
    if varyinhib==0
        c0(1) = c01*conc_fac(n);
        sub_conc(n) = log10(c0(1));
        c0(2) = c02;
    else
        c0(1) = c01;
        c0(2) = c02*conc_fac(n);
        inh_conc(n) = log10(c0(2));
    end
    m = 0;
    sol = pdepe(m,@us1pde,@us1ic,@us1bc,x,t);
    us = sol(:,:,1); % Extract first component as substrate
    ui = sol(:,:,2); % Extract second component as inhibitor
    usc = sol(:,:,3); % Extract third component as intracellular substrate
    uic = sol(:,:,4); % Extract fourth component as intracellular inhibitor
    for i = 1:ttt
        cell_sub_content(i) = trapz(x,usc(i,:))/Cellfac;
        cell_inh_content(i) = trapz(x,uic(i,:))/Cellfac;
    end
    % check conservation of substrate by initial and final amounts
    sub_init = trapz(x,us(1,:)) + trapz(x,usc(1,:))/Cellfac;
    sub_final = trapz(x,us(ttt,:)) + trapz(x,usc(ttt,:))/Cellfac;
    inh_init = trapz(x,ui(1,:)) + trapz(x,uic(1,:))/Cellfac;
    inh_final = trapz(x,ui(ttt,:)) + trapz(x,uic(ttt,:))/Cellfac;
    if sub_init > 0
        sub_error = (sub_final - sub_init)/sub_init;
        if abs(sub_error) > 0.001
            disp('*** Warning: conservation error of substrate');
        end
    end
    if inh_init > 0
        inh_error = (inh_final - inh_init)/inh_init;
        if abs(inh_error) > 0.001
            disp('*** Warning: conservation error of substrate');
        end
    end
    sub_uptake_rate(n) = 60*(cell_sub_content(ttt)-cell_sub_content(1))/t(ttt);
    inh_uptake_rate(n) = 60*(cell_inh_content(ttt)-cell_inh_content(1))/t(ttt);
end;
%***** end loop over substrate or inhibitor concentrations ******

figure(1);
title('Average rate of uptake: substrate and inhibitor');
if varyinhib==0
    plot(sub_conc,sub_uptake_rate,'-o',sub_conc,inh_uptake_rate,'-x')
    xlabel('log[Substrate conc. in microM]');
else
    plot(inh_conc,sub_uptake_rate,'-o',inh_conc,inh_uptake_rate,'-x')
    xlabel('log[Inhibitor conc. in microM]');
end
ylabel('Average uptake rate (nmol/cm2/min)');

fid = fopen('Uptake_vs_conc.txt','w');
fprintf(fid, 'Time = %f s\n', Tmax);
fprintf(fid, 'Initial concentration (uM) Average uptake nmol/cm2/min\n');
fprintf(fid, ' substrate   inhibitor    substrate   inhibitor\n');
for n = 1:nconcs
    fprintf(fid, '%12.8f %12.8f %12.8e %12.8e\n', c0(1),c0(2),...
        sub_uptake_rate,inh_uptake_rate);
end;

% **** plots and data files only for last set of parameters! ****

% concentrations at cell membrane (x=0) and in reservoir (x=L)
figure(2);
plot(t,us(:,1),'-o',t,us(:,length(x)),'-o',t,ui(:,1),'-x',t,ui(:,length(x)),'-x')
title('Substrate and inhibitor concentration at membrane, in reservoir');
xlabel('Time (s)');
ylabel('Concentration (microM)');
fid = fopen('Concentrations-time.txt','w');
fprintf(fid, '           Substrate (uM)                         Inhibitor (uM)\n');
fprintf(fid,...
'Time (s)   Membrane    Reservoir     Cells        Membrane    Reservoir   Cells\n');
for n = 1:length(t)
    fprintf(fid, '%f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n',...
    t(n),us(n,1),us(n,length(x)),usc(n,1),ui(n,1),ui(n,length(x)),uic(n,1));
end;

% intracellular accumulations of substrate and inhibitor
figure(3);
plot(t,cell_sub_content,'-o',t,cell_inh_content,'-x')
title('Cell content: substrate and inhibitor');
xlabel('Time (s)');
ylabel('Amount (nmol/cm2)');
fid = fopen('Uptake-time.txt','w');
fprintf(fid, 'Time (s) Substrate Inhibitor nmol/cm2\n');
for n = 1:length(t)
    fprintf(fid, '%f %12.8f %12.8f\n',...
        t(n),cell_sub_content(n),cell_inh_content(n));
end;

% % surface plot of substrate
% figure(4);
% surf(x,t,us);    
% title('Substrate concentration');
% xlabel('Distance (cm)');
% ylabel('Time (s)');
% 
% % surface plot of inhibitor
% figure(5);
% surf(x,t,ui);    
% title('Inhibitor concentration');
% xlabel('Distance (cm)');
% ylabel('Time (s)');
% 
% % surface plot of substrate
% figure(6);
% surf(x,t,usc);    
% title('Substrate concentration');
% xlabel('Distance (cm)');
% ylabel('Time (s)');
% 
% % surface plot of inhibitor
% figure(7);
% surf(x,t,uic);    
% title('Inhibitor concentration');
% xlabel('Distance (cm)');
% ylabel('Time (s)');
% % concentration profiles at various time points
% figure(8);
% plot(x,usc(1,:),x,usc(9,:),x,usc(14,:),x,usc(19,:),x,usc(22,:),...
%     x,usc(ttt,:))
% % axis([0 0.2 0 .02]);
% title('Inhibitor concentration profiles (t=0,3,6,12,30,60 s)');
% xlabel('Position (cm)');
% ylabel('Inhibitor concentration (microM)');
% % solution profile
% figure;
% hold all;
% for it = 1:length(t)
%     plot(x,u(it,:))
% end
% hold off;
% title('Concentration profiles');
% xlabel('Distance (cm)');
% ylabel('Concentration');

fclose('all');
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

function u0 = us1ic(x)
% initial solute level, solute level outside layer
% includes initial solute-free, ramped or smooth initial concentration 
global c0 L00 ramp
if x >= L00
  u0 = c0';
else
    if ramp == 0
        u0 = [0 0 0 0]';
    else
        u0 = (1 - (1 - x/L00)^ramp)*c0';
    end
end
end
    
% --------------------------------------------------------------------------

function [pl,ql,pr,qr] = us1bc(~,ul,~,~,~)
% Michaelis-Menten at x = 0, with inhibition for substrate
% zero flux at x = L
global Jm Kt Cellfac Jl1
Jsol = Jm(1)*ul(1)/(Kt(1)*(1 + ul(2)/Kt(2)) + ul(1))... 
    - Jm(3)*ul(3)/(Kt(3)*(1 + ul(4)/Kt(4)) + ul(3))...
    + Jl1*ul(1);    % include linear uptake term

Jinh = Jm(2)*ul(2)/(Kt(2)*(1 + ul(1)/Kt(1)) + ul(2)) ...
    - Jm(4)*ul(4)/(Kt(4)*(1 + ul(3)/Kt(3)) + ul(4));

pl = [Jsol Jinh -Jsol*Cellfac -Jinh*Cellfac]';
ql = [-1 -1 -1 -1]';
pr = [0 0 0 0]';
qr = [1 1 1 1]';
end
