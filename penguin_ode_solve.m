function [M,Cflux,Pflux] = penguin_ode_solve(Pe,beta,alpha,eps,N,steps,tstep,toly,problem,shape)
% Penguin ODE Solver
%   M = penguin_ode_solve(Pe,beta,alpha,eps,N,steps,tstep,toly,problem,shape)
%   produces the cell array M of size 1x(steps+1) of the free boundary evolution. 
%
%   [M, Cflux, Pflux] = penguin_ode_solve(Pe,beta,alpha,eps,N,steps,tstep,toly,problem,shape)
%   produces the cell array M, Cflux, Pflux for the free boundary evolution
%   and the canonical and physical heat fluxes of the final free boundary
%   shape, respectively.
%
%
% INPUTS
%   Pe          = Peclet number - strength of the flow.
%   beta        = thermal conductivity ratio - strength of interior effect.
%   alpha       = Poisson equation forcing.
%   eps         = amplitude of noise.
%   N           = series truncation of conformal map. 
%                   Number of unknowns = n = 2N+3.
%   steps       = number of time steps.
%   tstep       = value of each time step.
%   toly        = relative and absolute tolerances on ode15i.
%   problem     = problem choice. The ode solver is applicable to six free
%                   boundary problems depending on the inclusion/exclusion
%                   of the exterior, interior and area conserving terms in
%                   the Stefan boundary condition. See [1] for more
%                   details.
%                    0 (ext); 1 (int); 2 (ext+int); 3 (ext+areaC), 
%                       4 (int+areaC), 5 (ext+int+areaC).
%   shape       = initial shape of the free boundary. 
%                   See penguin_initial_shape for all shape choices.
%
% OUTPUTS
%   M           = cell array of size 1x(steps+1) of free boundary data.
%                   Each entry (t-1) gives a list of complex z=x+iy values 
%                   for points with coordinates (x,y) on the free boundary 
%                   at time step t. (t=0: initial shape.)
%
% OPTIONAL OUTPUTS
%   Cflux       = heat flux in the canonical plane across the (upper half 
%                   of the) free boundary of the final shape. A 1x2 cell 
%                   array; Cflux{1} is a list of theta values of the boundary 
%                   points in the canonical plane; Cflux{2} is the canonical 
%                   heat flux at each point.
%
%   Pflux       = heat flux in the physical plane across the (upper half 
%                   of the) free boundary of the final shape. A 1x2 cell 
%                   array; Pflux{1} is a list of z values of the boundary 
%                   points in the physical plane; Cflux{2} is the physical 
%                   heat flux at each point.
%
% REFERENCES
%   [1]	        Harris, S.J., McDonald, N.R. (2023) "Penguin Huddling: 
%               A Continuum Model". Acta Appl. Math. 185, 7. 
%               https://doi.org/10.1007/s10440-023-00578-2.
%
% END OF DOCUMENTATION
%
%Code
warning('off','MATLAB:rankDeficientMatrix'); warning('off','CHEBFUN:aaa:Froissart'); warning('off','MATLAB:polyshape:repairedBySimplify'); %unnecessary warnings turned off

%Fixed Parameters
n=2*N+3; opts = odeset('RelTol',toly, 'AbsTol',toly);  %n = number of unknown coefficients = 2N+3; opts = tolerances for ode15i,
tmin=0; tmax = tstep*steps; tsteps=steps+1; T = linspace(tmin,tmax,tsteps); %tsteps = total number of steps (including t=0); T=time vector,
[c0_est, cdot0_est] = penguin_initial_shape(shape,N,n); %estimates for initial coefficients (c0) and their time derivative (cdot0)
randFun = -1+2*rand(1,n+1);  %random (spatial) function: assigns each position theta a value between -1 and 1

%f, decic and ode15i
f = @(t,c,cdot)odeqns(t,c,cdot,N,n,Pe,beta,alpha,eps,randFun,problem); tic; %Defining ftn f based on set of ODEs -- see odeqns ftn
[c0,cdot0] = decic(f,tmin,c0_est,1,cdot0_est,[],opts); %initial conditions for ode15i - c0_est is entirely fixed, cdot0_est is entirely non-fixed
[~,coeff] = ode15i(f,T,c0,cdot0); coeff=coeff.'; %ode solver

%Creating the polygons of the free boundary
Pts=1000; Theta=linspace(0,2*pi,Pts); Zeta = exp(1i*Theta); %Pts= number of points on the polygon
for j=1:tsteps
z=coeff(1,j)*Zeta;
    for kk=2:N+1
    z=z+(coeff(kk,j)+1i*coeff(N+1+kk,j))*Zeta.^(-(kk-2)); %At step j, complex Laurent coefficient c_(kk-2) = coeff(kk) + 1i*coeff(kk+(N+1)) 
    end
M{j}=z; %output as a cell array M
end

%Finding the heat flux of final shape - canonical and physical planes
Theta2 = linspace(0,2*pi,n+1); Zeta2=exp(1i*Theta2); %Use a smaller number of points
z2=coeff(1,end)*Zeta2; dzeta2=coeff(1,end); %z, dz/dzeta
    for kk=2:N+1
    z2=z2+(coeff(kk,end)+1i*coeff(N+1+kk,end))*Zeta2.^(-(kk-2));
    dzeta2=dzeta2-(kk-2).*(coeff(kk,end)+1i*coeff(N+1+kk))*Zeta2.^(-(kk-2)-1);
    end
zedeze2=Zeta2.*dzeta2; %zeta*dz/dzeta
centr = coeff(2,end)+1i*coeff(N+3,end); %conformal centre of polygon

fsigma=sigsolve(Theta2,Pe,n,eps,randFun); %sigma of final shape
fFd =(AAA_LS_solve(z2.',N,alpha,centr)).';
fomega= real(zedeze2.*(fFd-alpha*conj(z2)/2)); %omega of final shape

cflux = fsigma-beta*fomega; Cflux = {Theta2, cflux}; %theta and canonical heat flux
pflux = (cflux)./abs(dzeta2); Pflux = {centrepoly(z2), pflux}; %z and physical heat flux

totaltime = toc; fprintf('ODE Solver Complete. \n Total runtime: %.1f seconds.\n',totaltime); %displays runtime
end
