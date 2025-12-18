%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%   "Modelling collective migration of phenotypically heterogeneous   %%%
%%%          cell populations: from single-cell dynamics                %%%     
%%%                to population-level behaviours"                      %%%
%%%                                                                     %%%
%%%        T. Lorenzi, N. Loy (*), L. Preziosi, C. Villa, 2026          %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%  Code for the numerical integration of the microscopic (2.1)-(2.2)  %%%
%%%  with a Monte Carlo scheme in 1D.  [copyright: Nadia Loy (*)]       %%%
%%%                                                                     %%%
%%%  This version differs from the one in 'MICRO_SimMC_1D_tot.m'        %%%
%%%  because it also takes into account the interactions in which       %%%
%%%  both phenotypic switching and directional changes occur - an       %%%
%%%  effect of order Dt^2.                                              %%%
%%%                                                                     %%%
%%% (*) nadia.loy@polito.it                                             %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all 
close all
clc


%% Parameters

%%% Parameters from general set up
par = Parameters();
a = par.x1min;
b = par.x1max;
c = par.x2min;
d = par.x2max;

%%% Choose scenario (A) or (B) for ymaxL
% par.ymaxL = par.ymaxLA; % A) 1
par.ymaxL = par.ymaxLB; % B) 0

%%% Choose Kappa definition:
Kappa = 'DD'; % Dirac delta
% Kappa = 'VM'; % Von Mises

%%% Physical domain
Dx1 = .01;
Dx2 = .01;
xg1 = [a+0.5*Dx1:Dx1:b-0.5*Dx1];
xg2 = [c+0.5*Dx2:Dx2:d-0.5*Dx2];
Nx1 = length(xg1);
Nx2 = length(xg2);
xg1 = [xg1,b+0.5*Dx1];
xg2 = [xg2,d+0.5*Dx2];

%%% Velocity domain
Vmax = 1.7;
dv = 1e-2;
vmod = [0:dv:Vmax];

%%% Phenotype domain
ymax = 1;
ymin = 0;
dy = 0.05;
yy = ymin:dy:ymax;
Ny = length(yy);

%%% Time discretization
Dt = 1e-4;
Tmax = 12;
t = 0:Dt:Tmax;
Nt = length(t);

%%% Frequencies
mu = 1;     % rate of velocity change
lambda = 1; % rate of phenotypic change 

%%% Monte Carlo parameter
N = 1e6;    % Number of particles


%% Set up and initial conditions

%%% ECM density (directional cue)
Mat=@(s) par.Mmin+par.Mgr*(s-c);

%%% Macroscopic cell density rho_0 (as for macro simulation)
rho_0 = zeros(Nx2,1);
L2 = 1;
rho0=@(s2) max(0,1*(s2<-L2 & s2 >= 2*L2-par.x1max)+(s2>=-L2).*(1-((s2+L2)./L2).^2)+(s2<2*L2-par.x1max).*(1-((2*L2-par.x1max-s2)./L2).^2));
rho_0(xg2<-L2) = 1;
% rhoA(1:Nx1,x2<par.x2min_cells) = 0;
rho_0(xg2>=-L2) = 1-(((xg2(xg2>=-L2)+L2)./L2).^2);
rho_0(xg2<2*L2-par.x1max) = 1-(((2*L2-par.x1max-xg2(xg2<2*L2-par.x1max))./L2).^2);

rho_0(rho_0<0) = 0;
mass=sum(rho_0)*Dx2;

%%% Sampling for the Monte Carlo simulation from rho_0
Nr=N;
i0nr=[1:N];
while Nr>1
    r20=c+(d-c)*rand(Nr,1);
    y120=max(max(rho0(xg2(1:end-1))))*rand(Nr,1);
    inr=find(y120<rho0(r20));
    x20(i0nr(inr))=r20(inr);
    i0nr(inr)=[];
    Nr=length(i0nr);
end
x2=x20'; % Initial positions

%%% Initial velocity 
vx20=zeros(N,1);
vx2=vx20;

%%% Initial phenotype
y0=zeros(N,1); 
y=y0;

%% Time iterations

it=1;
for nt=1:Nt

    %%% Random permutation of the agents
    rn=randperm(N);
    x2=x2(rn);
    vx2=vx2(rn);
    y=y(rn);
    
    %%% Interacting agents of the current itaration
    Tv=binornd(1,(mu/par.eps)*Dt,N,1);
    Tl=binornd(1,(lambda/par.eps)*Dt,N,1);
    Ttot=Tv+Tl;
    itot=find(Ttot==2);
    iv=find(Tv);
    il=find(Tl);
    Tv((ismember(iv,itot)))=0;
    Tl((ismember(il,itot)))=0;
    iv=find(Tv);
    il=find(Tl);

    %%% Direction dynamics: Rejection method
    nr=length(iv);
    ivnr=iv;
    ic=0;
    pr_dir=zeros(N,1);
    i=1;
    while nr>1
        i=i+1;
        r1=binornd(1,0.5,nr,1);
        r1(r1==0)=-1;
        r2=(Mat(min(d*ones(nr,1),max((c)*ones(nr,1),x2(ivnr)+par.Rmax)))./(2*(0.1+x2(ivnr)-c))).*rand(nr,1);
        iin=find(r2<M(min(d*ones(nr,1),max((c)*ones(nr,1),x2(ivnr)+(par.Rmin+(par.Rmax-par.Rmin).*y(ivnr)).*r1)))./(2*(0.1+x2(ivnr)-c)));
        pr_dir(ivnr(iin))=r1(iin);
        
        ivnr(iin)=[];
        nr=length(ivnr);
    end
    
    dir2=pr_dir;
    vx2(iv)=dir2(iv)*par.vmaxL;

    %%% Phenotypic switch
    y(il)=par.ymaxL;

    %%% Agents doing both things
    % Phenotypic switch
    y(itot)=par.ymaxL;
    % Direction dynamics: Rejection method
    nr=length(itot);
    ivnr=itot;
    ic=0;
    pr_dir=zeros(N,1);
    i=1;
    while nr>1
        i=i+1;
        r1=binornd(1,0.5,nr,1);
        r1(r1==0)=-1;
        r2=(Mat(min(d*ones(nr,1),max((c)*ones(nr,1),x2(ivnr)+par.Rmax)))./(2*(0.1+x2(ivnr)-c))).*rand(nr,1);
        iin=find(r2<Mat(min(d*ones(nr,1),max((c)*ones(nr,1),x2(ivnr)+(par.Rmin+(par.Rmax-par.Rmin).*y(ivnr)).*r1)))./(2*(0.1+x2(ivnr)-c)));
        pr_dir(ivnr(iin))=r1(iin);
        
        ivnr(iin)=[];
        nr=length(ivnr);
    end

    dir2=pr_dir;
    vx2(itot)=dir2(itot)*par.vmaxL;
   
    %%%Transport 
    x2=x2+Dt*vx2;

    %%%'Boundary conditions'
    iy_c=find(x2<c+Dx2);
    iy_d=find(x2>d-Dx2);
    
    vx2(find(vx2(iy_c)<0))=-vx2(find(vx2(iy_c)<0));
    vx2(find(vx2(iy_d)>0))=-vx2(find(vx2(iy_d)>0));

    x2(iy_c)=c+(Dx2).*rand(length(iy_c),1);
    x2(iy_d)=d-Dx2+(Dx2).*rand(length(iy_d),1);

    %%% Check they are all inside the domain
    check=(x2>c) & (x2<d);
    
    %%% Save at 12h and 24h
    if Dt*nt==12
        x2c12=x2;
    elseif Dt*nt==24
        x2c24=x2;
    end
    
    disp(sprintf('Time=%f,Mass=%f',Dt*nt,length(find(check))))
    
end

%% Save setup and results
folder = 'Saved_Data/'; 
save([folder,'datasMC_1D_tot_eps',num2str(par.eps),'.mat'])

