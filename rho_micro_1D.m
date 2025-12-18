%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%   "Phenotype-structuring of non-local kinetic models of cell        %%%     
%%%           migration driven by environmental sensing"                %%%
%%%                                                                     %%%
%%%              T. Lorenzi, N. Loy, C. Villa, 2026                     %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%  Code to reconstruct cell density from microscopic model results    %%%
%%%                in 1D  [copyright: Nadia Loy (*)]                    %%%
%%%                                                                     %%%
%%% (*) nadia.loy@polito.it                                             %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear xg2
clear rho
clear rhoc0
clear rhoc12
clear rhoc24

Dx2=0.15;

xg2=[c+0.5*Dx2:Dx2:d-0.5*Dx2];
Nx2=length(xg2);
xg2=[xg2,d+0.5*Dx2];

x20=x20';

for j=1:Nx2
    den= (x20>=xg2(j)-0.5*Dx2) & (x20<=xg2(j+1)-0.5*Dx2);
    nj=length(find(den));
        if den==zeros(length(den),1)
            nj=0;
        end
    rhoc0(j)=nj/(N*Dx2);
end    

rhoc0=rhoc0/sum(rhoc0*Dx2);

rhoc0=mass*rhoc0;


for j=1:Nx2
    den=(x2c12>=xg2(j)-0.5*Dx2) & (x2c12<=xg2(j+1)-0.5*Dx2);
    nj=length(find(den));
    if den==zeros(length(den),1)
        nj=0;
    end
    rhoc12(j)=nj/(N*Dx2); 
end
rhoc12=rhoc12/sum(rhoc12*Dx2);
rhoc12=mass*rhoc12;
%rhoc24=(LM)*abs(-20+0.5)*rhoc24;

for j=1:Nx2
    den=(x2c24>=xg2(j)-0.5*Dx2) & (x2c24<=xg2(j+1)-0.5*Dx2);
    nj=length(find(den));
    if den==zeros(length(den),1)
        nj=0;
    end
    rhoc24(j)=nj/(N*Dx2); 
end
rhoc24=rhoc24/sum(rhoc24*Dx2);
rhoc24=mass*rhoc24;
