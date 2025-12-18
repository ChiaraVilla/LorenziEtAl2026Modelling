%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%     Approximation of velocity term in nonlocal 1D sensing region    %%%
%%%             [copyright: Chiara Villa (*), Nadia Loy]                %%%
%%%                                                                     %%%
%%% (*) chiara[dot]villa[at]math[dot]cnrs[dot]fr                        %%%
%%%                                                                     %%%
%%%   Disclaimer: this is not optimised (i.e. overall code will be slow %%%
%%%   if this needs to be re-calculated at every iteration in time),    %%%
%%%   for improvements see Gerisch 2010 (doi:10.1093/imanum/drp027)     %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [Ue2A,Ue2B,De2A,De2B] = Nonlocal_advection_1D_DD(par)
%%%
%%% UT given at cell edges (not cell centers)
%%%
%%% For now: code given chosing gamma = Dirac delta cenetered in R(y)
%%% Kappa given by Dirac delta distribution

a = par.x1min;
b = par.x1max;
c = par.x2min;
d = par.x2max;

% x1 and x2 discretisation at cell edges (maintaining same dx)
x1A = 0.5*(par.str1+a); % Middle of LM stripe
x1B = 0.5*(par.str2+par.str1); %Middle of FN stripe
x2 = (par.x2min):par.dx2:(par.x2max); 
Nx1 = 1; % 1D
Nx2 = length(x2);

% y discretisation (choose dy)
dy = 0.1;
y = (0:dy:1);
Ny = length(y);

% vhat discretisation
vhat = [-1,1];
th = [pi/2,3*pi/2];
dth = 1;

% ECM density gradient (determining direction of motion) - def (SM5.3)
Mat=@(s2) par.Mmin+par.Mgr*(s2-c);    
% Constant of integration (only if gamma = Dirac delta) > x2-dependent
cM = repmat(1./((2)*(par.Mmin+par.Mgr*(x2-c))),Nx1,1);  

% Mean post-reorientation speed - i.e. u_psi in def (4.3)
vbar_f=@(s1,s2) par.vmaxL*((s1<par.str1) | (s1 >par.str2)) + par.vmaxF*((s1>=par.str1) & (s1 <=par.str2));

% Mean post-transition phenotype - def (4.7)
% For the purpose of the 1D test: choose ybarLB for both and just evaluate
% A at LN stripe and B at FN stripe
ybar_fA=@(s1,s2) par.ymaxLB*((s1<par.str1) | (s1 >par.str2)) + par.ymaxF*((s1>=par.str1) & (s1 <=par.str2));
ybar_fB=@(s1,s2) par.ymaxLB*((s1<par.str1) | (s1 >par.str2)) + par.ymaxF*((s1>=par.str1) & (s1 <=par.str2));

% Sensing radius - def (4/6)
R_f=@(yd) par.Rmin + yd.*(par.Rmax-par.Rmin); 

% Initialise UT at cell edges
UTx2A = zeros(Nx1,Nx2);
UTx2B = zeros(Nx1,Nx2);
DTx2A = zeros(Nx1,Nx2);
DTx2B = zeros(Nx1,Nx2);

% Evaluate ybar outside the loop 
ybA = repmat(ybar_fA(x1A,x2),1,Nx2);
ybB = repmat(ybar_fB(x1B,x2),1,Nx2);

for j=1:2

    % Evaluate vbA at each (x1,x2) + R(ybar)*(cos(th),sin(th))
    % > avoid exiting the domain
    vbA=vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1A+R_f(ybA).*cos(th(j)))),...
        min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybA).*sin(th(j))))); 
    vbB=vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1B+R_f(ybB).*cos(th(j)))),...
        min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybB).*sin(th(j))))); 
    
    % Add integration contribution to UT
    UTx2A=UTx2A+cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybA).*sin(th(j)))))*sin(th(j)).*vbA*dth;
    UTx2B=UTx2B+cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybB).*sin(th(j)))))*sin(th(j)).*vbB*dth;

end

for j=1:2

    % Evaluate vbA at each (x1,x2) + R(ybar)*(cos(th),sin(th))
    % > avoid exiting the domain
    vbA=vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1A+R_f(ybA).*cos(th(j)))),...
        min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybA).*sin(th(j))))); 
    vbB=vbar_f(min(b*ones(Nx1,Nx2),max(a*ones(Nx1,Nx2),x1B+R_f(ybB).*cos(th(j)))),...
        min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybB).*sin(th(j))))); 
        
    % Add integration contribution to D_T 
    DTx2A=DTx2A+cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybA).*sin(th(j))))).*((vbA.*sin(th(j))-UTx2B).^2)*dth;
    DTx2B=DTx2B+cM.*Mat(min(d*ones(Nx1,Nx2),max(c*ones(Nx1,Nx2),x2+R_f(ybB).*sin(th(j))))).*((vbB.*sin(th(j))-UTx2B).^2)*dth;

end


%%% Ghost points to compute \div UT (ensuring no-flux boundary conditions)
UTx2gpA = [UTx2A(:,1),UTx2A,UTx2A(:,end)];
UTx2gpB = [UTx2B(:,1),UTx2B,UTx2B(:,end)];
DTx2gpA = [DTx2A(:,1),DTx2A,DTx2A(:,end)];
DTx2gpB = [DTx2B(:,1),DTx2B,DTx2B(:,end)];

%%% UT_eps = UT( 1 - eps \div UT) 
Ue2A = UTx2A.*( 1 - par.eps*((UTx2gpA(:,3:end)-UTx2gpA(:,1:end-2))./(2*par.dx2)));
Ue2B = UTx2B.*( 1 - par.eps*((UTx2gpB(:,3:end)-UTx2gpB(:,1:end-2))./(2*par.dx2)));

%%% DT_eps = eps * DT
De2A = par.eps*DTx2A;
De2B = par.eps*DTx2B;


end               