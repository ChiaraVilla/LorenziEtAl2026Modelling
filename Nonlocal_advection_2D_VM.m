%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%   "Phenotype-structuring of non-local kinetic models of cell        %%%     
%%%           migration driven by environmental sensing"                %%%
%%%                                                                     %%%
%%%              T. Lorenzi, N. Loy, C. Villa, 2026                     %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%     Approximation of velocity term in nonlocal 2D sensing region    %%%

%%%             [copyright: Chiara Villa (*), Nadia Loy (**)]           %%%
%%%                                                                     %%%
%%% (*) chiara[dot]villa[at]math[dot]cnrs[dot]fr                        %%%
%%% (**) nadia.loy@polito.it                                            %%%
%%%                                                                     %%%
%%%   Disclaimer: this is not optimised (i.e. overall code will be slow %%%
%%%   if this needs to be re-calculated at every iteration in time),    %%%
%%%   for improvements see Gerisch 2010 (doi:10.1093/imanum/drp027)     %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ue1A,Ue2A,Ue1B,Ue2B,De11A,De12A,De22A,De11B,De12B,De22B] = Nonlocal_advection_2D_VM(par)
%%%
%%% UT_eps and DT_eps given at cell edges (corners)
%%%
%%% Code given chosing: 
%%% - gamma(r) = Dirac delta centered in R(y) > definition (2.12)
%%% - psi(v) = Dirac delta centered in u_\psi > satisfies (2.13)
%%% - M(x) = Mmin + Mgr( x2 - Lm ) > definition (SM5.3)
%%%
%%% Kappa(y): Von Mises (VM) distribution

k = par.ky; % Von Mises coefficient

a = par.x1min;
b = par.x1max;
c = par.x2min;
d = par.x2max;

dx1 = par.dx1;
dx2 = par.dx2;

% x1 and x2 discretisation at cell edges (maintaining same dx)
x1 = (par.x1min):par.dx1:(par.x1max);
x2 = (par.x2min):par.dx2:(par.x2max); 
Nx1 = length(x1);
Nx2 = length(x2);

% y discretisation (choose dy)
dy = 0.1;
y = (0:dy:1);
Ny = length(y);

% theta discretisation (choose Nth to ensure equally spaced points)
Nth = 60;
dth = 2*pi/Nth;
th = (0:dth:2*pi-dth)';

% ECM density gradient (determining direction of motion) 
Mat=@(s2) par.Mmin+par.Mgr*(s2-c);    
% Constant of integration (only if gamma = Dirac delta) > x2-dependent   
cM = repmat(1./((2*pi)*(par.Mmin+par.Mgr*(x2-c))),Nx1,1);    

% Mean post-reorientation speed - i.e. u_psi in def (4.3)
vbar_f=@(s1,s2) par.vmaxL*((s1<par.str1) | (s1 >par.str2)) + par.vmaxF*((s1>=par.str1) & (s1 <=par.str2));

% Mean post-transition phenotype - def (4.7)
ybar_fA=@(s1,s2) par.ymaxLA*((s1<par.str1) | (s1 >par.str2)) + par.ymaxF*((s1>=par.str1) & (s1 <=par.str2));
ybar_fB=@(s1,s2) par.ymaxLB*((s1<par.str1) | (s1 >par.str2)) + par.ymaxF*((s1>=par.str1) & (s1 <=par.str2));

% Sensing radius - def (4.6)
R_f=@(yd) par.Rmin + yd.*(par.Rmax-par.Rmin); 

% Evaluate ybar at each (x1,x2)
ybA=repmat(ybar_fA(x1',x2),1,Nx2);
ybB=repmat(ybar_fB(x1',x2),1,Nx2);

% Initialise uT, aT, sT, UT, AT, CT at cell edges
uTx1A = zeros(Nx1,Nx2);
uTx2A = uTx1A;
aT11A = uTx1A;
aT12A = uTx1A;
aT22A = uTx1A;
sTx1A = uTx1A;
sTx2A = uTx1A;
UTx1A = uTx1A;
UTx2A = uTx1A;
AT11A = uTx1A;
AT12A = uTx1A;
AT22A = uTx1A;
CT11A = uTx1A;
CT12A = uTx1A;
CT22A = uTx1A;

uTx1B = zeros(Nx1,Nx2);
uTx2B = uTx1B;
aT11B = uTx1B;
aT12B = uTx1B;
aT22B = uTx1B;
sTx1B = uTx1B;
sTx2B = uTx1B;
UTx1B = uTx1B;
UTx2B = uTx1B;
AT11B = uTx1B;
AT12B = uTx1B;
AT22B = uTx1B;
CT11B = uTx1B;
CT12B = uTx1B;
CT22B = uTx1B;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:Ny
    for j=1:Nth

        % Evaluate vbA at each (x1,x2) + R(y)*(cos(th),sin(th))
        % > avoid exiting the domain
        vbA = vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1'+R_f(y(i)).*cos(th(j)))),...
            min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(y(i)).*sin(th(j)))));

        vbB = vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1'+R_f(y(i)).*cos(th(j)))),...
            min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(y(i)).*sin(th(j)))));

        % Add integration contribution to uT
        MatijA = cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
            x2+R_f(y(i)).*sin(th(j)))));
        MatijB = cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),...
            x2+R_f(y(i)).*sin(th(j)))));

        uTx1A = uTx1A + MatijA*cos(th(j)).*vbA*dth;
        uTx2A = uTx2A + MatijA*sin(th(j)).*vbA*dth;

        uTx1B = uTx1B + MatijB*cos(th(j)).*vbB*dth;
        uTx2B = uTx2B + MatijB*sin(th(j)).*vbB*dth;

        % Add integration contribution to aT
        aT11A = aT11A + MatijA*(cos(th(j)).^2).*(vbA.^2)*dth;
        aT12A = aT12A + MatijA*(cos(th(j)).*sin(th(j))).*(vbA.^2)*dth;
        aT22A = aT22A + MatijA*(sin(th(j)).^2).*(vbA.^2)*dth;

        aT11B = aT11B + MatijB*(cos(th(j)).^2).*(vbB.^2)*dth;
        aT12B = aT12B + MatijB*(cos(th(j)).*sin(th(j))).*(vbB.^2)*dth;
        aT22B = aT22B + MatijB*(sin(th(j)).^2).*(vbB.^2)*dth;
    end

    % Von Mises distribution
    VMiA = (0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybA)...
        + 2*pi-(pi-2*pi.*ybA).*(y(i)-2*pi.*ybA)))*dy;
    VMiB = (0.5/besseli(0,k)).*exp(-k*cos(-(pi-2*pi.*ybB)...
        + 2*pi-(pi-2*pi.*ybB).*(y(i)-2*pi.*ybB)))*dy;

    % Evaluate UT integrating uT against Von Mises distribution
    UTx1A = UTx1A + uTx1A.*VMiA;
    UTx2A = UTx2A + uTx2A.*VMiA;

    UTx1B = UTx1B + uTx1B.*VMiB;
    UTx2B = UTx2B + uTx2B.*VMiB;

    % Evaluate AT integrating aT against Von Mises distribution
    AT11A = AT11A + aT11A.*VMiA;
    AT12A = AT12A + aT12A.*VMiA;
    AT22A = AT22A + aT22A.*VMiA;

    AT11B = AT11B + aT11B.*VMiB;
    AT12B = AT12B + aT12B.*VMiB;
    AT22B = AT22B + aT22B.*VMiB;

    % Evaluate C integrating uT⊗uT against Von Mises distribution
    CT11A = CT11A + (uTx1A.^2).*VMiA;
    CT12A = CT12A + (uTx1A.*uTx2A).*VMiA;
    CT22A = CT22A + (uTx2A.^2).*VMiA;

    CT11B = CT11B + (uTx1B.^2).*VMiB;
    CT12B = CT12B + (uTx1B.*uTx2B).*VMiB;
    CT22B = CT22B + (uTx2B.^2).*VMiB;

    % Evaluate c integrating \div(uT)uT against Von Mises distribution
    sTx1A = sTx1A + div(uTx1A,uTx2A,dx1,dx2).*uTx1A.*VMiA;
    sTx2A = sTx2A + div(uTx1A,uTx2A,dx1,dx2).*uTx2A.*VMiA;

    sTx1B = sTx1B + div(uTx1B,uTx2B,dx1,dx2).*uTx1B.*VMiB;
    sTx2B = sTx2B + div(uTx1B,uTx2B,dx1,dx2).*uTx2B.*VMiB;

end

% Compute variance-covariace matrix DT = AT - UT⊗UT
DT11A = AT11A - UTx1A.^2;
DT12A = AT12A - UTx1A.*UTx2A;
DT22A = AT22A - UTx2A.^2;

DT11B = AT11A - UTx1A.^2;
DT12B = AT12A - UTx1A.*UTx2A;
DT22B = AT22A - UTx2A.^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% UT_eps = UT(1- eps/2 \div UT) - eps/2 cT - \eps \div DT
Ue1A = UTx1A.*( 1 - 0.5*par.eps*div(UTx1A,UTx2A,dx1,dx2) ) ...
    - 0.5*par.eps*sTx1A - par.eps*div(DT11A,DT12A,dx1,dx2);
Ue2A = UTx2A.*( 1 - 0.5*par.eps*div(UTx1A,UTx2A,dx1,dx2) ) ...
    - 0.5*par.eps*sTx2A - par.eps*div(DT12A,DT22A,dx1,dx2);

Ue1B = UTx1B.*( 1 - 0.5*par.eps*div(UTx1B,UTx2B,dx1,dx2) ) ...
    - 0.5*par.eps*sTx1A - par.eps*div(DT11A,DT12A,dx1,dx2);
Ue2B = UTx2B.*( 1 - 0.5*par.eps*div(UTx1B,UTx2B,dx1,dx2) ) ...
    - 0.5*par.eps*sTx2B - par.eps*div(DT12B,DT22B,dx1,dx2);

%%% DT_eps = \eps DT +eps/2 (CT - UT ⊗ UT )
De11A = par.eps*DT11A + 0.5*par.eps*(CT11A - UTx1A.^2);
De12A = par.eps*DT12A + 0.5*par.eps*(CT12A - UTx1A.*UTx2A);
De22A = par.eps*DT22A + 0.5*par.eps*(CT22A - UTx2A.^2);

De11B = par.eps*DT11B + 0.5*par.eps*(CT11B - UTx1B.^2);
De12B = par.eps*DT12B + 0.5*par.eps*(CT12B - UTx1B.*UTx2B);
De22B = par.eps*DT22B + 0.5*par.eps*(CT22B - UTx2B.^2);

end               

function divu = div(u1,u2,dx1,dx2)
    divu = 0*u1;
    divu(2:end-1,2:end-1) = ( u1(3:end,2:end-1) - u1(1:end-2,2:end-1) )./(2*dx1) ...
        + ( u2(2:end-1,3:end) - u2(2:end-1,1:end-2) )./(2*dx2);
end