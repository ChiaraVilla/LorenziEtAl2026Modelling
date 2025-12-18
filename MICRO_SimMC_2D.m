%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%   "Phenotype-structuring of non-local kinetic models of cell        %%%     
%%%           migration driven by environmental sensing"                %%%
%%%                                                                     %%%
%%%              T. Lorenzi, N. Loy, C. Villa, 2026                     %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%  Code for the numerical integration of the microscopic (2.1)-(2.2)  %%%
%%%  with a Monte Carlo scheme in 2D.  [copyright: Nadia Loy (*)]       %%%
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

%%% Average velocity imposed by the stripes
vbar=@(s1,s2) par.vmaxL*((s1<par.str1) | (s1>par.str2))+par.vmaxF*((s1>=par.str1) & (s1<=par.str2));

%%% Average phenotype imposed by the stripes
ybar=@(s1,s2) par.ymaxL*((s1<par.str1) | (s1>par.str2))+par.ymaxF*((s1>=par.str1) & (s1<=par.str2));

%%% ECM density (directional cue)
Mat=@(s) par.Mmin+par.Mgr*(s-c);    

%%% Macroscopic cell density rho_0 (as for macro simulation)
rho_0 = zeros(Nx1,Nx2);
L2 = 1;
rho0=@(s2) max(0,1*(s2<-L2 & s2 >= 2*L2-par.x1max)+(s2>=-L2).*(1-((s2+L2)./L2).^2)+(s2<2*L2-par.x1max).*(1-((2*L2-par.x1max-s2)./L2).^2));
rho_0(1:Nx1,xg2<-L2) = 1;

for i=1:Nx1
    rho_0(i,xg2>=-L2) = 1-(((xg2(xg2>=-L2)+L2)./L2).^2);
    rho_0(i,xg2<2*L2-par.x1max) = 1-(((2*L2-par.x1max-xg2(xg2<2*L2-par.x1max))./L2).^2);
end
rho_0(rho_0<0) = 0;
mass = sum(sum(rho_0)*Dx1*Dx2);

%%% Sampling for the Monte Carlo simulation from rho_0
x10=(b)*rand(N,1);
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
x1=x10;  % Initial positions 
x2=x20'; % Initial positions

%%% Initial velocity 
vx10=zeros(N,1);
vx20=zeros(N,1);
vx1=vx10;
vx2=vx20;

%%% Initial phenotype 
y0=zeros(N,1); 
y=y0;

%% Time iterations

it=1;
for nt=1:Nt

    %%% Random permutation of the agents
    rn=randperm(N);
    x1=x1(rn);
    x2=x2(rn);
    vx1=vx1(rn);
    vx2=vx2(rn);
    y=y(rn);
    
    %%% Interacting agents of the current iteration
    Tv=binornd(1,(mu/par.eps)*Dt,N,1);
    Tl=binornd(1,(lambda/par.eps)*Dt,N,1);
    iv=find(Tv);
    il=find(Tl);
    
    %%% Phenotype switch
    switch Kappa

        case 'DD' % Dirac delta
            y(il)=ybar(x1(il),x2(il));

        case 'VM' % Von Mises
            %%%%% DOWNLOAD the function 'vmrand.m' by Dylan Muir (2024). vmrand(fMu, fKappa, varargin) 
            %%%%% (https://www.mathworks.com/matlabcentral/fileexchange/37241-vmrand-fmu-fkappa-varargin), 
            %%%%% MATLAB Central File Exchange. Retrieved October 22, 2024.
            yb=ybar(x1(il),x2(il));
            alpha=vmrand((2*pi*yb)-pi,10,length(iv),1);
            y(il)=(alpha+pi)./(2*pi);
    end
    
    %%% Direction switch: Rejection method
    nr=length(iv);
    ivnr=iv;
    ic=0;
    pr_dir=zeros(N,1);
    i=1;
    while nr>1
        i=i+1;
        r1=(2*pi)*rand(nr,1);
        r2=(Mat(min(d*ones(nr,1),max((c)*ones(nr,1),x2(ivnr)+par.Rmax)))./(2*pi*(0.1+x2(ivnr)-c))).*rand(nr,1);
        iin=find(r2<Mat(min(d*ones(nr,1),max((c)*ones(nr,1),x2(ivnr)+(par.Rmin+(par.Rmax-par.Rmin).*y(ivnr)).*sin(r1))))./(2*pi*(0.1+x2(ivnr)-c)));
        pr_dir(ivnr(iin))=r1(iin);
        ivnr(iin)=[];
        nr=length(ivnr);
    end
   
    dir1=cos(pr_dir(iv));
    dir2=sin(pr_dir(iv));
    
    vx1(iv)=dir1.*(vbar(min(max(a*ones(length(iv),1),x1(iv)+dir1.*(par.Rmin+(par.Rmax-par.Rmin).*y(iv))),b*ones(length(iv),1)),...
        min(max(c*ones(length(iv),1),x2(iv)+dir2.*(par.Rmin+(par.Rmax-par.Rmin).*y(iv))),d*ones(length(iv),1))));
    vx2(iv)=dir2.*(vbar(min(max(a*ones(length(iv),1),x1(iv)+dir1.*(par.Rmin+(par.Rmax-par.Rmin).*y(iv))),b*ones(length(iv),1)),...
        min(max(c*ones(length(iv),1),x2(iv)+dir2.*(par.Rmin+(par.Rmax-par.Rmin).*y(iv))),d*ones(length(iv),1))));
   
    %%% Transport
    x1=x1+Dt*vx1;
    x2=x2+Dt*vx2;

    %%% Boundary conditions
    ix_a=find(x1<a+Dx1);
    ix_b=find(x1>b-Dx1);
    iy_c=find(x2<c+Dx2);
    iy_d=find(x2>d-Dx2);
    
    vx1(find(vx1(ix_a)<0))=-vx1(find(vx1(ix_a)<0));
    vx1(find(vx1(ix_b)>0))=-vx1(find(vx1(ix_b)>0));
    vx2(find(vx2(iy_c)<0))=-vx2(find(vx2(iy_c)<0));
    vx2(find(vx2(iy_d)>0))=-vx2(find(vx2(iy_d)>0));

    x1(ix_a)=a+(Dx1).*rand(length(ix_a),1);
    x1(ix_b)=b-Dx1+(Dx1).*rand(length(ix_b),1);
    x2(iy_c)=c+(Dx2).*rand(length(iy_c),1);
    x2(iy_d)=d-Dx2+(Dx2).*rand(length(iy_d),1);

    %%% Check they are all inside the domain
    check=(x1>a) & (x1<b) & (x2>c) & (x2<d);
    
    %%% Save at 24h and 48h
    if Dt*nt==24
        x1c24=x1;
        x2c24=x2;
    elseif Dt*nt==48
        x1c48=x1;
        x2c48=x2;
    end
    
    disp(sprintf('Time=%f,Mass=%f',Dt*nt,length(find(check))))
    
end

%% Save setup and results
folder = 'Saved_Data/'; 
save([folder,'datasMC_2D_eps',num2str(par.eps),'_N',num2str(N),'.mat'])

