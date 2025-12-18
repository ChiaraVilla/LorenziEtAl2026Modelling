%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%               MUSCL** scheme for computing the flux                 %%%
%%%                  [copyright: Chiara Villa (*)]                      %%%
%%%                                                                     %%%
%%% (*) chiara[dot]villa[at]math[dot]cnrs[dot]fr                        %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%%  ** Monotonic Upwind-Centered Scheme for Conservation Laws          %%%
%%%     - Flux conserving                                               %%%
%%%     - Piecewise linear reconstruction                               %%%
%%%     - Total variation diminishing (TVD)                             %%%
%%%     - Uniform grid                                                  %%%
%%%     - Explicit (flux calculated at known timestep t_n)              %%%
%%%       or with Correction term (flux calculated at t_{n+1/2})        %%%
%%%       in case the drift velocity is NOT a function of time          %%%
%%%     - Account for more than one dimension                           %%%
%%%     - Impose No-flux Boundary Conditions                            %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dfdx] = MUSCL_GP(u,v,nx,ic1,dx,dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Equation: du/dt + df/dx = 0, f=uv
%%%
%%% INPUT:
%%% - u: advected quantity at cell centers (+ 4 ghost points)
%%% - v: drift velocity at cell centers or interfaces
%%% - nx: number of grid cells (excluding ghost points)
%%% - ic1: index cell 1 in u vector (to account for ghost points)
%%% - dx: grid cell width (assuming uniform grid)
%%% - dt: timestep width (optional: only needed for correction term)
%%%
%%% OUTPUT:
%%% - dfdx: df/dx computed at cell centers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Check u has the correct size: nx + at least 2 ghost points on each side
if size(u,ndims(u))<nx+4
    msg = 'Incorrect size for advected quantity array: MUSCL requires at least 2 ghost points at each end of the domain';
    error(msg)
end

%%% Approximate v at cell interfaces i-1/2 (from 1/2 to nx+1/2), size nx+1
switch size(v,ndims(v))
    case size(u,ndims(u)) % given at cell centers, including all ghost points
        vim12 = 0.5*(v(:,ic1-1:ic1+nx-1)+v(:,ic1:ic1+nx));
    case nx+1 % given at cell interfaces of main grid       
        if ndims(v)==2 && size(v,1)==size(u,1)+1
            vim12 = 0.5*(v(1:end-1,:)+v(2:end,:));
        else
            vim12 = v;
        end
    case nx+2 % given at cell centers, including two ghost points (l+r)
        vim12 = 0.5*(v(:,1:end-1)+v(:,2:end));
    case size(u,ndims(u))+1 % given at cell interfaces of whole grid, including ghost points
        vim12 = v(:,ic1:ic1+nx);
    otherwise
        msg = 'Incorrect size for velocity array: check MUSCL function requirements';
        error(msg)
end


%%%% Ratio of consecutive gradients at i-1/2 interface (from 1/2 to nx+3/2), size nx+2
rim12 = (u(:,ic1-1:ic1+nx)-u(:,ic1-2:ic1+nx-1))./(u(:,ic1:ic1+nx+1)-u(:,ic1-1:ic1+nx));
rim12(isnan(rim12)) = 0; % avoid NaN values where u is constant

%%% Correction term: if dv/dt=0 and we want df/dx at t_n+1/2
switch nargin
    case 6
        cim12 = abs(vim12.*(dt/dx)); 
        %cim12 = 0; % Explicit scheme (df/dx at t_n) has no correction term
    case 5
        cim12 = 0; % No correction term without dt given as input
    otherwise
        error('Not enough input arguments in MUSCL_GP function')
end

%%% Approximate u at i-1/2 from left and right, using a flux limiter, size nx+1
uLim12 = u(:,ic1-1:ic1+nx-1) + 0.5*(1-cim12).*fl(rim12(:,1:nx+1)).*(u(:,ic1:ic1+nx)-u(:,ic1-1:ic1+nx-1)); 
uRim12 = u(:,ic1:ic1+nx) - 0.5*(1-cim12).*fl(rim12(:,2:nx+2)).*(u(:,ic1+1:ic1+nx+1)-u(:,ic1:ic1+nx));

%%% Flux at i-1/2 (from 1/2 to nx+1/2)
fim12 = uLim12.*max(0,vim12) + uRim12.*min(0,vim12); 

%%% No-flux Boundary Conditions (at 1/2 and nx+1/2)
fim12(:,1) = 0;
fim12(:,nx+1) = 0;

%%% df/dx at cell centers (1st order centered difference) 
dfdx = (fim12(:,2:nx+1)-fim12(:,1:nx))./dx; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Flux limiters 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [flr] = fl(r)

    %%% Donor cell / Upwind
    flr = 0;
    
    %%% Lax-Wendroff
    % flr = 1; 
    
    %%% minimod 
    % flr = max(0,min(1,r));
    
    %%% ospre
    % flr = 1.5.*((r).^2+r)./(r.^2+r+1);
    % flr(isfinite(flr)~=1) = 1.5;
    
    %%% Koren
    % flr = max(0,min(2*r,min((1+2*r)./3,2)));
    
    %%% superbee
    % flr = max(0,max(min(2*r,1),min(r,2)));
    
    %%% monotonized central (MC)
    % flr = max(0,min(2*r,min(0.5*(1+r),2)));
    
    %%% van Leer
    % flr = (r+abs(r))./(1+abs(r));
    % flr(isfinite(r)~=1) = 2;

end