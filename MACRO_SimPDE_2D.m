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
%%%       over laminin & fibronectin, with the macroscopic model of     %%%
%%%           Eq.(3.48)  [copyright: Chiara Villa (*)]                  %%%
%%%                                                                     %%%
%%% (*) chiara[dot]villa[at]math[dot]cnrs[dot]fr                        %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear all
close all
set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 10)
set(0,'defaultaxeslinewidth',1)
set(0,'defaultpatchlinewidth',1)
set(0,'defaultlinelinewidth',2)
set(0,'defaultTextInterpreter','latex')
%%% To save data
date = '111225'; 
folder = 'Saved_Data/'; 
add = ''; % default

%% Parameters and discretisation

%%% Choose whether to use Kappa = Dirac delta (DD) or von Mises (VM)
%Kappa = 'DD';
Kappa = 'VM';

%%% Problem parameters (in units of 10^-2 cm)
par = Parameters();
par.eps = 0.001;
add = ['_eps',num2str(par.eps)]; % any additional specification for data file titles

%%% Space discretisation
dx = 0.05; % Choose grid refinement
par.dx1 = dx; 
x1 = (par.x1min+0.5*par.dx1):par.dx1:(par.x1max-0.5*par.dx1); % Cell centres
Nx1 = length(x1);
par.dx2 = dx; 
x2 = (par.x2min+0.5*par.dx2):par.dx2:(par.x2max-0.5*par.dx2); % Cell centres
Nx2 = length(x2);

%%% Final time
Tend = 48;

%% Initial conditions
rho0 = zeros(Nx1,Nx2);
L2 = 1;
rho0(1:Nx1,x2<-L2) = 1;
for i=1:Nx1
    rho0(i,x2>=-L2) = 1-(((x2(x2>=-L2)+L2)./L2).^2);
    rho0(i,x2<2*L2-par.x1max) = 1-(((2*L2-par.x1max-x2(x2<2*L2-par.x1max))./L2).^2);
end
rho0(rho0<0) = 0;
% Calculate initial mass to check mass conservation
MassA = sum(rho0,"all")*par.dx2*par.dx1; 
MassB = MassA;

%% Computation of the advection velocity
%  > Details of the problem-specific definitions affecting UT and DT are 
%    given in  the file 'Nonlocal_advection.m' 
%  > UT and DT are calculated and given at the cell edges (corners)
%  > UT_eps = UT(1- eps/2 \div UT) - eps/2 cT - \eps \div DT
%  > DT_eps = \eps DT +eps/2 (CT - UT ⊗ UT )

switch Kappa
    case 'DD'
        [UTx1A,UTx2A,UTx1B,UTx2B,DT11A,DT12A,DT22A,DT11B,DT12B,DT22B] = Nonlocal_advection_2D_DD(par);
    case 'VM'
        [UTx1A,UTx2A,UTx1B,UTx2B,DT11A,DT12A,DT22A,DT11B,DT12B,DT22B] = Nonlocal_advection_2D_VM(par);
    otherwise
        error('Choose appropriate definition for Kappa (DD or VM)')
end

% save([folder,'Saved_',date,'_Setup_',Kappa,add,'.mat'],'par','Kappa','UTx2A',"UTx1A","UTx1B","UTx2B",'DT11A','DT12A','DT22A','DT11B','DT12B','DT22B')

%% Setup solver for equation (3.48) 
%%% Equation:
%%% \dt p + \div ( p * UT ) = \div  ( DT_eps * \div p )
%%% UT_eps = UT(1- eps/2 \div UT) - eps/2 cT - \eps \div DT
%%% DT_eps = \eps DT +eps/2 (CT - UT ⊗ UT )

%%% UT_eps components at Cell Edges Horizontal/Vertical 
UTx1cehA = 0.5*( UTx1A(:,1:end-1)+UTx1A(:,2:end) );
UTx2cevA = 0.5*( UTx2A(1:end-1,:)+UTx2A(2:end,:) );

UTx1cehB = 0.5*( UTx1B(:,1:end-1)+UTx1B(:,2:end) );
UTx2cevB = 0.5*( UTx2B(1:end-1,:)+UTx2B(2:end,:) );

%%% DT_eps components at Cell Edges Horizontal/Vertical 
DT11cehA = 0.5*( DT11A(:,2:end)+DT11A(:,1:end-1) );
DT12cehA = 0.5*( DT12A(:,2:end)+DT12A(:,1:end-1) );
DT21cevA = 0.5*( DT12A(2:end,:)+DT12A(1:end-1,:) );
DT22cevA = 0.5*( DT22A(2:end,:)+DT22A(1:end-1,:) );

DT11cehB = 0.5*( DT11B(:,2:end)+DT11B(:,1:end-1) );
DT12cehB = 0.5*( DT12B(:,2:end)+DT12B(:,1:end-1) );
DT21cevB = 0.5*( DT12B(2:end,:)+DT12B(1:end-1,:) );
DT22cevB = 0.5*( DT22B(2:end,:)+DT22B(1:end-1,:) );

%%% Explicit gradient matrices
GMx1 = GM_def(Nx1,par.dx1);
GMx2 = GM_def(Nx2,par.dx2);


%% Iterate in time (MUSCL scheme)

%%% Time discretisation & the CFL condition
dt = 0.01; % Choose maximum dt
UmaxA = max(max(max(UTx1A)),max(max(UTx2A)));
UmaxB = max(max(max(UTx1B)),max(max(UTx2B)));
Umax = max(UmaxA,UmaxB);
DmaxA = max(max(max(DT11A)),max(max(DT22A)));
DmaxB = max(max(max(DT11B)),max(max(DT22A)));
Dmax = max(DmaxA,DmaxB);
CFL_adv = dx/Umax;
CFL_diff = (dx^2)/(4*Dmax);
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
rhostoreA = zeros((Nt-1)*dt,size(rho0,1),size(rho0,2));
rhostoreB = rhostoreA;
MassA = zeros(1,Nt);
MassB = zeros(1,Nt);
rhoA = rho0;
rhoB = rho0;

%%% Compute recurrent factors for computing of ghost points 
% % consistently with no-flux boundary conditions
% % NOTE: the contribution of DT12/DT21 is excluded because really small
gp1x1A = (DT11cehA(1,:)./par.dx1 - UTx1cehA(1,:)./2)./(DT11cehA(1,:)./par.dx1 + UTx1cehA(1,:)./2);
gpNx1A = (DT11cehA(end,:)./par.dx1 + UTx1cehA(end,:)./2)./(DT11cehA(end,:)./par.dx1 - UTx1cehA(end,:)./2);
gp1x2A = (DT22cevA(:,1)./par.dx2 - UTx2cevA(:,1)./2)./(DT22cevA(:,1)./par.dx2 + UTx2cevA(:,1)./2);
gpNx2A = (DT22cevA(:,end)./par.dx2 + UTx2cevA(:,end)./2)./(DT22cevA(:,end)./par.dx2 - UTx2cevA(:,end)./2);

gp1x1B = (DT11cehB(1,:)./par.dx1 - UTx1cehB(1,:)./2)./(DT11cehB(1,:)./par.dx1 + UTx1cehB(1,:)./2);
gpNx1B = (DT11cehB(end,:)./par.dx1 + UTx1cehB(end,:)./2)./(DT11cehB(end,:)./par.dx1 - UTx1cehB(end,:)./2);
gp1x2B = (DT22cevB(:,1)./par.dx2 - UTx2cevB(:,1)./2)./(DT22cevB(:,1)./par.dx2 + UTx2cevB(:,1)./2);
gpNx2B = (DT22cevB(:,end)./par.dx2 + UTx2cevB(:,end)./2)./(DT22cevB(:,end)./par.dx2 - UTx2cevB(:,end)./2);


%%% Time iteration
for i=1:Nt-1

    %%% ADVECTION
    %%% Compute \div(f) given f=UT*rho, using MUSCL: add ghost points
    rhogp1A = [zeros(1,size(rhoA,2)); rhoA(1,:).*gp1x1A  ;rhoA;...
        rhoA(end,:).*gpNx1A  ;zeros(1,size(rhoA,2))];
    dfdx1A = MUSCL_GP(rhogp1A',UTx1cehA',Nx1,3,par.dx1,dt)';
    rhogp2A = [zeros(size(rhoA,1),1),rhoA(:,1).*gp1x2A,...
        rhoA, rhoA(:,end).*gpNx2A, zeros(size(rhoA,1),1)];
    dfdx2A = MUSCL_GP(rhogp2A,UTx2cevA,Nx2,3,par.dx2,dt);

    rhogp1B = [zeros(1,size(rhoB,2)); rhoB(1,:).*gp1x1B  ;rhoB;...
        rhoB(end,:).*gpNx1B  ;zeros(1,size(rhoB,2))];
    dfdx1B = MUSCL_GP(rhogp1B',UTx1cehB',Nx1,3,par.dx1,dt)';
    rhogp2B = [zeros(size(rhoB,1),1),rhoB(:,1).*gp1x2B,...
        rhoB, rhoB(:,end).*gpNx2B, zeros(size(rhoB,1),1)];
    dfdx2B = MUSCL_GP(rhogp2B,UTx2cevB,Nx2,3,par.dx2,dt);

    %%% ANISOTROPIC DIFFUSION
    %%% Compute div(S) with S=DT*\grad(rho) 
    S1cehA = DT11cehA.*Gcehx1(rhoA,par.dx1,Nx1,Nx2) + ...
        DT12cehA.*Gcehx2(rhoA,par.dx2,Nx1,Nx2); % Only needed at Horizontal CE
    S2cevA = DT21cevA.*(Gcehx2(rhoA',par.dx1,Nx2,Nx1)') + ...
        DT22cevA.*(Gcehx1(rhoA',par.dx2,Nx2,Nx1)'); % Only needed at Vertical CE
    dSdxA = GMx1*S1cehA + (GMx2*(S2cevA'))';

    S1cehB = DT11cehB.*Gcehx1(rhoB,par.dx1,Nx1,Nx2) + ...
        DT12cehB.*Gcehx2(rhoB,par.dx2,Nx1,Nx2); % Only needed at Horizontal CE
    S2cevB = DT21cevB.*(Gcehx2(rhoB',par.dx1,Nx2,Nx1)') + ...
        DT22cevB.*(Gcehx1(rhoB',par.dx2,Nx2,Nx1)'); % Only needed at Vertical CE
    dSdxB = GMx1*S1cehB + (GMx2*(S2cevB'))';

    %%% PDE
    rhoA = rhoA - dt*(dfdx1A+dfdx2A + dSdxA );

    rhoB = rhoB - dt*(dfdx1B+dfdx2B + dSdxB );

    %%% Store Mass to check Mass Conserrvation
    MassA(i) = sum(rhoA,"all")*par.dx2*par.dx1;
    MassB(i) = sum(rhoB,"all")*par.dx2*par.dx1;

    %%% Store every hour
    if mod(i-1,1/dt)==0
        rhostoreA(1+(i-1)*dt,:,:) = rhoA;
        rhostoreA(1+(i-1)*dt,:,:) = rhoB;
        subplot(1,2,1)
        Plot_StripeS_Timet(rhoA,x1,x2,par,par.x1min,par.x1max,par.x2min,par.x2max,max(max(rhoA)),'n')
        title(['$t=$',num2str(t(i))],'Interpreter','latex')
        subplot(1,2,2)
        Plot_StripeS_Timet(rhoB,x1,x2,par,par.x1min,par.x1max,par.x2min,par.x2max,max(max(rhoB)),'n')
        title(['$t=$',num2str(t(i))],'Interpreter','latex')
        drawnow
    end

end

% save([folder,'Saved_',date,'_Results_',Kappa,add,'.mat'],'rhostoreA','rhostoreB','MassA','MassB','rhoA','rhoB','rho0')

%% Final plot
% Kappa = 'VM';
% load([folder,'Saved_',date,'_Setup_',Kappa,add,'.mat'])
% load([folder,'Saved_',date,'_Results_',Kappa,add,'.mat'])

%%% Plot Stripes
figure(1)
Plot_stripes(rho0,rhoA,rhoB,par,x1,x2)

% %%% Plot to check mass conservation
% figure(2)
% subplot(1,2,1)
% plot(MassA)
% title('Mass $\rho_A$','Interpreter','latex')
% subplot(1,2,2)
% plot(MassB)
% title('Mass $\rho_B$','Interpreter','latex')

%% Gradient matrices (finite volume scheme with no-flux BCs)

function GM = GM_def(Nx,dx)
%%% Simple 1D Gradient matrix (d_x) :
%%% - multiplies quantity given at cell edges in direction of gradient
%%% - product returs directional gradient at cell centers

    %%% Coefficient matrix
    GM = zeros(Nx,Nx+1);
    GM(1:Nx,1:Nx) = GM(1:Nx,1:Nx) + eye(Nx);
    GM(1:Nx,2:Nx+1) = GM(1:Nx,2:Nx+1) - eye(Nx);
    GM = GM./dx;

end

function Gcehx1 = Gcehx1(p,dx1,Nx1,Nx2) 
%%% Function to compute the gradient in x1 direction of p:
%%% - input p given at cell centers
%%% - output d_x1(p) given at Horizontal Cell Edges  | . |

    %%% Initialize matrix 
    Gcehx1 = zeros(Nx1+1,Nx2);
    
    %%% Internal entries 
    Gcehx1(2:Nx1,:) = p(2:Nx1,:)-p(1:Nx1-1,:);
    %%% Boundary terms (due to ghost point ensuring no-flux BCs are satisfied)
    Gcehx1(1,:) = zeros(1,Nx2);
    Gcehx1(Nx1+1,:) = zeros(1,Nx2);
    
    %%% Output
    Gcehx1 = Gcehx1./dx1;

end

function Gcehx2 = Gcehx2(p,dx2,Nx1,Nx2) 
%%% Function to compute the gradient in x2 direction of p:
%%% - input p given at cell centers
%%% - output d_x2(p) given at Horizontal Cell Edges  | . |

    %%% Extend p with ghost points
    pgp = zeros(Nx1+2,Nx2+2);
    pgp(2:Nx1+1,2:Nx2+1) = p;
    pgp(1,:) = pgp(2,:);
    pgp(Nx1+2,:) = pgp(Nx1+1,:);
    pgp(:,1) = pgp(:,2);
    pgp(:,Nx2+2) = pgp(:,Nx2+1);
    
    %%% Entries using ghost points
    Gcehx2 = 0.25*( pgp(1:Nx1+1,1:Nx2) + pgp(1:Nx1+1,3:Nx2+2) + ...
                    pgp(2:Nx1+2,1:Nx2) + pgp(2:Nx1+2,3:Nx2+2) );
    %%% Output
    Gcehx2 = Gcehx2./dx2;

end

