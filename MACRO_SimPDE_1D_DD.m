%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%   "Phenotype-structuring of non-local kinetic models of cell        %%%     
%%%           migration driven by environmental sensing"                %%%
%%%                                                                     %%%
%%%              T. Lorenzi, N. Loy, C. Villa, 2026                     %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%  Code to simulate Stripe migration essays of Goodman et al. (1989)  %%%
%%%      over laminin & fibronectin, with the macroscopic model of      %%%
%%%  Eq.(3.50) in 1D for K = Dirac delta [copyright: Chiara Villa (*)]  %%%
%%%                                                                     %%%
%%% (*) chiara[dot]villa[at]math[dot]cnrs[dot]fr                        %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 14)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultTextInterpreter','latex')
%%% To save data
date = '110924'; 
folder = 'Saved_Data/'; 
add = ''; % default

%% Parameters 

%%% Problem parameters (in units of 10^-2 cm)
par = Parameters();
par.eps = 0.001;
D = num2str(par.eps);

%%% Space discretisation
dx = 0.01;                   % Choose grid refinement
Nx1 = 1;
par.dx2 = dx; 
x2 = (par.x2min+0.5*par.dx2):par.dx2:(par.x2max-0.5*par.dx2); % Cell centres
Nx2 = length(x2);

%%% Final time
Tend = 40;

%% Initial conditions
rho0 = zeros(1,Nx2);
L2 = 1;
rho0(1,x2<-L2) = 1;
rho0(1,x2>=-L2) = 1-(((x2(x2>=-L2)+L2)./L2).^2);
rho0(1,x2<2*L2-par.x2max) = 1-(((2*L2-par.x2max-x2(x2<2*L2-par.x2max))./L2).^2);
rho0(rho0<0) = 0;
% Calculate initial mass to check mass conservation
Mass = sum(rho0,"all")*par.dx2; 

%% Computation of the advection velocity
%  > Details of the problem-specific definitions affecting UT are given in 
%    the file 'Nonlocal_advection.m' 
%  > UT_eps and DT_eps are calculated and given at the cell edges
%  > UT_eps = UT(1-eps\div UT)
%  > DT_eps = eps*DT

% 1D: A = LN stripe and B = FN stripe
[UTx2A,UTx2B,DTx2A,DTx2B] = Nonlocal_advection_1D_DD(par);

% save([folder,'Saved_',date,'_Setup_1D_eps',D,'.mat'],'par','UTx2A',"UTx2B",'DTx2A',"DTx2B")

%% Setup solver for equation (3.50) 
%%%
%%% \dt p + \div ( p * UT_eps ) = \div \div ( DT_eps * p )

%%% DT_eps needed at cell centers in this formulation of the equation
DTA_c = 0.5*(DTx2A(2:end)+DTx2A(1:end-1));
DTB_c = 0.5*(DTx2B(2:end)+DTx2B(1:end-1));

%%% Explicit diffusion matrix: \div \div (...)
DMx2 = DM_def(Nx2,par.dx2);

%% Iterate in time (MUSCL scheme)

%%% CFL condition
dt = 0.01; % Choose maximum dt
Umax = max(max(max(UTx2A)),max(max(UTx2B)));
Dmax = max(max(max(DTx2A)),max(max(DTx2B)));
CFL_adv = dx/Umax;
CFL_diff = (dx^2)/(2*Dmax);
if dt>CFL_adv || dt>CFL_diff
   dt = 0.5*min(CFL_adv,CFL_diff);
   % CFL condition changed the dt, ensure solutions will be stored hourly
   for i=floor(1.0/dt):-1:1
       if mod(1,1.0/i)==0
           dt = 1.0/i;
           break;
       end
   end
end
t = 0:dt:Tend;
Nt = length(t);

%%% Initialise storing arrays and dynamics variables
rhostoreA = zeros((Nt-1)*dt,length(rho0));
rhostoreB = rhostoreA;
MassA = zeros(1,Nt);
MassB = zeros(1,Nt);
rhoA = rho0;
rhoB = rho0;

%%% Compute recurrent factors for computing of ghost points 
% % consistently with no-flux boundary conditions
gp1x2A = (DTx2A(1)./par.dx2 - UTx2A(1)./2)./(DTx2A(1)./par.dx2 + UTx2A(1)./2);
gpNx2A = (DTx2A(end)./par.dx2 - UTx2A(end)./2)./(DTx2A(end)./par.dx2 + UTx2A(end)./2);

gp1x2B = (DTx2B(1)./par.dx2 - UTx2B(1)./2)./(DTx2B(1)./par.dx2 + UTx2B(1)./2);
gpNx2B = (DTx2B(end)./par.dx2 - UTx2B(end)./2)./(DTx2B(end)./par.dx2 + UTx2B(end)./2);


%%% Iterate

for i=1:Nt

    %%% Compute df/dx given f=UT*rho, using MUSCL: add ghost points
    rhogp2A = [0,gp1x2A*rhoA(1),rhoA,gpNx2A*rhoA(end),0];
    dfdx2A = MUSCL_GP(rhogp2A,UTx2A,Nx2,3,par.dx2,dt);

    rhogp2B = [0,gp1x2B*rhoB(1),rhoB,gpNx2B*rhoB(end),0];
    dfdx2B = MUSCL_GP(rhogp2B,UTx2B,Nx2,3,par.dx2,dt);

    %%% Explicit in time
    rhoA = rhoA - dt*(dfdx2A) + dt*(DMx2*(DTA_c.*(rhoA))')';

    rhoB = rhoB - dt*(dfdx2B) + dt*(DMx2*(DTB_c.*(rhoB))')';


    %%% Store Mass to check Mass Conserrvation
    MassA(i) = sum(rhoA,"all")*par.dx2;
    MassB(i) = sum(rhoB,"all")*par.dx2;

    %%% Store every hour
    if mod(i-1,1/dt)==0
        rhostoreA(round(1+(i-1)*dt),:) = rhoA;
        rhostoreA(round(1+(i-1)*dt),:) = rhoB;
        subplot(2,1,1)
        plot(rhoA )
        title('rho (LN)')
        subplot(2,1,2)
        plot(rhoB )
        title('rho (FN)')
        drawnow
     end

end

% save([folder,'Saved_',date,'_Results_1D_eps',D,'.mat'],'rhostoreA','rhostoreB','Mass','rhoA','rhoB','rho0')

%% Final plot
% load(['Saved_',date,'_Setup.mat'])
% load(['Saved_',date,'_Results.mat'])

%%% Plot stripes
figure(1)
subplot(1,2,1)
plot(x2,rho0,'k:')
hold on
plot(x2,rhoA,'r')
title('LN')
subplot(1,2,2)
plot(x2,rho0,'k:')
hold on
plot(x2,rhoB,'b')
title('FN')

%%% Plot to check mass conservation
figure(2)
plot(MassA)

%% Diffusion matrix (finite volume scheme with no-flux BCs)

function DM = DM_def(Nx,dx)

    %%% Coefficient matrix
    DM = -2*eye(Nx);
    DM(2:Nx,1:Nx-1) = DM(2:Nx,1:Nx-1) + eye(Nx-1);
    DM(1:Nx-1,2:Nx) = DM(1:Nx-1,2:Nx) + eye(Nx-1);

    %%% No-flux boundary conditions (from finite volume scheme)
    DM(1,1) = -1;
    DM(1,2) = 1;
    DM(Nx,Nx-1) = 1;
    DM(Nx,Nx) = -1;

    %%% Output
    DM = DM./(dx^2);

end