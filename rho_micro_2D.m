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
%%%                in 2D  [copyright: Nadia Loy (*)]                    %%%
%%%                                                                     %%%
%%% (*) nadia.loy@polito.it                                             %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear xg1
clear xg2
clear rho
clear rhoc0
clear rhoc24
clear rhoc48

xg1=[a+0.5*Dx1:Dx1:b-0.5*Dx1];
xg2=[c+0.5*Dx2:Dx2:d-0.5*Dx2];
Nx1=length(xg1);
Nx2=length(xg2);
xg1=[xg1,b+0.5*Dx1];
xg2=[xg2,d+0.5*Dx2];

x20=x20';

for i=1:Nx1
    for j=1:Nx2
        den=(x10>=xg1(i)-0.5*Dx1) & (x10<xg1(i+1)-0.5*Dx1) & (x20>=xg2(j)-0.5*Dx2) & (x20<=xg2(j+1)-0.5*Dx2);
        nij=length(find(den));
            if den==zeros(length(den),1)
                nij=0;
            end
        rhoc0(i,j)=nij/(N*Dx1*Dx2);
    end    
end
rhoc0=rhoc0/sum(sum(rhoc0*Dx1*Dx2));
rhoc0=mass*rhoc0;

for i=1:Nx1
    for j=1:Nx2
        den=(x1c24>=xg1(i)-0.5*Dx1) & (x1c24<xg1(i+1)-0.5*Dx1) & (x2c24>=xg2(j)-0.5*Dx2) & (x2c24<=xg2(j+1)-0.5*Dx2);
        ni=length(find(den));
        if den==zeros(length(den),1)
            ni=0;
        end
        rhoc24(i,j)=ni/(N*Dx1*Dx2);
    end    
end
rhoc24=rhoc24/sum(sum(rhoc24*Dx1*Dx2));
rhoc24=mass*rhoc24;

for i=1:Nx1
    for j=1:Nx2
        den=(x1c48>=xg1(i)-0.5*Dx1) & (x1c48<xg1(i+1)-0.5*Dx1) & (x2c48>=xg2(j)-0.5*Dx2) & (x2c48<=xg2(j+1)-0.5*Dx2);
        ni=length(find(den));
        if den==zeros(length(den),1)
            ni=0;
        end
        rhoc48(i,j)=ni/(N*Dx1*Dx2);
    end    
end
rhoc48=rhoc48/sum(sum(rhoc48*Dx1*Dx2));
rhoc48=mass*rhoc48;
