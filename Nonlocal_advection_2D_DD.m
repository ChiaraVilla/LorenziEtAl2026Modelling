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
%%%                      for Kappa = Dirac delta                        %%%
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

function [Ue1A,Ue2A,Ue1B,Ue2B,De11A,De12A,De22A,De11B,De12B,De22B] = Nonlocal_advection_2D_DD(par)
%%%
%%% UT given at cell edges, DT given at cell centers
%%%
%%% Code given chosing: 
%%% - gamma(r) = Dirac delta centered in R(y) > definition (2.12)
%%% - psi(v) = Dirac delta centered in u_\psi > satisfies (2.13)
%%% - M(x) = Mmin + Mgr( x2 - Lm ) > definition (SM5.3)
%%%
%%% Kappa(y): Dirac delta (DD) distribution

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

% Initialise UT and DT at cell edges
UTx1A = zeros(Nx1,Nx2);
UTx2A = zeros(Nx1,Nx2);
DT11A = zeros(Nx1,Nx2);
DT12A = zeros(Nx1,Nx2);
DT22A = zeros(Nx1,Nx2);

UTx1B = zeros(Nx1,Nx2);
UTx2B = zeros(Nx1,Nx2);
DT11B = zeros(Nx1,Nx2);
DT12B = zeros(Nx1,Nx2);
DT22B = zeros(Nx1,Nx2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluate ybar at each (x1,x2)
ybA = repmat(ybar_fA(x1',x2),1,Nx2);
ybB = repmat(ybar_fB(x1',x2),1,Nx2);

% Initialise computation arrays (for speed)
vbAs = zeros(Nth,Nx1,Nx2);
vbBs = vbAs;
cMAs = vbAs;
cMBs = vbAs;

for j=1:Nth

    % Evaluate vbA at each (x1,x2) + R(ybar)*(cos(th),sin(th))
    % > avoid exiting the domain
    vbA = vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1'+R_f(ybA).*cos(th(j)))),...
        min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybA).*sin(th(j))))); 
    vbAs(j,1:Nx1,1:Nx2) = vbA;

    vbB = vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1'+R_f(ybB).*cos(th(j)))),...
        min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybB).*sin(th(j))))); 
    vbBs(j,1:Nx1,1:Nx2) = vbB;

    % Evaluate integrand contribution from ECM
    cMA = cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybA).*sin(th(j)))));
    cMAs(j,1:Nx1,1:Nx2) = cMA;

    cMB = cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybB).*sin(th(j)))));
    cMBs(j,1:Nx1,1:Nx2) = cMB;

    % Add integration contribution to UT
    UTx1A = UTx1A + cMA*cos(th(j)).*vbA*dth;
    UTx2A = UTx2A + cMA*sin(th(j)).*vbA*dth;

    UTx1B = UTx1B + cMB*cos(th(j)).*vbB*dth;
    UTx2B = UTx2B + cMB*sin(th(j)).*vbB*dth;

end

for j=1:Nth

    cMA(1:Nx1,1:Nx2) = cMAs(j,1:Nx1,1:Nx2);
    vbA(1:Nx1,1:Nx2) = vbAs(j,1:Nx1,1:Nx2);

    cMB(1:Nx1,1:Nx2) = cMBs(j,1:Nx1,1:Nx2);
    vbB(1:Nx1,1:Nx2) = vbBs(j,1:Nx1,1:Nx2);

    % Add integration contribution to DT
    DT11A = DT11A + cMA.*((vbA.*cos(th(j))-UTx1A).^2)*dth;
    DT12A = DT12A + cMA.*(vbA.*cos(th(j))-UTx1A).*(vbA.*sin(th(j))-UTx2A)*dth;
    DT22A = DT22A + cMA.*((vbA.*sin(th(j))-UTx2A).^2)*dth;

    DT11B = DT11B + cMB.*((vbB.*cos(th(j))-UTx1B).^2)*dth;
    DT12B = DT12B + cMB.*(vbB.*cos(th(j))-UTx1B).*(vbB.*sin(th(j))-UTx2B)*dth;
    DT22B = DT22B + cMB.*((vbB.*sin(th(j))-UTx2B).^2)*dth;

end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% UT_eps = UT(1- eps \div UT) - \eps \div DT
Ue1A = UTx1A.*( 1 - par.eps*div(UTx1A,UTx2A,dx1,dx2) ) ...
    - par.eps*div(DT11A,DT12A,dx1,dx2);;
Ue2A = UTx2A.*( 1 - par.eps*div(UTx1A,UTx2A,dx1,dx2) )...
    - par.eps*div(DT11A,DT12A,dx1,dx2);;

Ue1B = UTx1B.*( 1 - par.eps*div(UTx1B,UTx2B,dx1,dx2) )...
    - par.eps*div(DT11B,DT12B,dx1,dx2);;
Ue2B = UTx2B.*( 1 - par.eps*div(UTx1B,UTx2B,dx1,dx2) )...
    - par.eps*div(DT11B,DT12B,dx1,dx2);;

%%% DT_eps = \eps DT 
De11A = par.eps*DT11A;
De12A = par.eps*DT12A;
De22A = par.eps*DT22A;

De11B = par.eps*DT11B;
De12B = par.eps*DT12B;
De22B = par.eps*DT22B;

end

function divu = div(u1,u2,dx1,dx2)
    divu = 0*u1;
    divu(2:end-1,2:end-1) = ( u1(3:end,2:end-1) - u1(1:end-2,2:end-1) )./(2*dx1) ...
        + ( u2(2:end-1,3:end) - u2(2:end-1,1:end-2) )./(2*dx2);
end
