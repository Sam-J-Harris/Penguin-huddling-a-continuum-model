function [sigma,omega,C] = SOCcalc(theta,z,zedze,Pe,beta,alpha,eps,randFun,c,N,n,centr,problem)
% Calculating Sigma, Omega and C Function
%   [sigma,omega,C] = SOCcalc(theta,z,zedze,Pe,beta,alpha,eps,randFun,c,N,n,centr,problem)
%   finds the sigma, omega and C terms in the Polubarinova-Galin (PG) type 
%   equation. See [1] for more information.
%
% INPUTS
%   theta       = azimuthal angle of points around the zeta-disk.
%   z           = list of complex z=x+iy values representing points with 
%                   coordinates (x,y) on the free boundary.
%   zedze       = zeta*(dz/dzeta) -- see odeqns.
%   Pe          = Peclet number - strength of the flow.
%   beta        = thermal conductivity ratio - strength of interior effect.
%   alpha       = Poisson equation forcing.
%   eps         = amplitude of noise.
%   randFun     = random (spatial) function: assigns each position theta a 
%                   value between -1 and 1.
%   c           = initial Laurent coefficient.
%   N           = series truncation of conformal map. 
%   n           = number of unknowns = 2N+3.
%   centr       = (conformal) centre of the shape traced by the free boundary.
%   problem     = problem choice. The ode solver is applicable to six free
%                   boundary problems depending on the inclusion/exclusion
%                   of the exterior, interior and area conserving terms in
%                   the Stefan boundary condition. See [1] for more
%                   details.
%                    0 (ext); 1 (int); 2 (ext+int); 3 (ext+areaC), 
%                       4 (int+areaC), 5 (ext+int+areaC).
%
% OUTPUTS
%   sigma       = exterior effect term in the PG equation.
%   omega       = interior effect term in the PG equation.
%   C           = area conserving term in the PG equation.
%
% REFERENCES
%   [1]	        Harris, S.J., McDonald, N.R. (2023) "Penguin Huddling: 
%               A Continuum Model". Acta Appl. Math. 185, 7. 
%               https://doi.org/10.1007/s10440-023-00578-2.
%
% END OF DOCUMENTATION
%
%Code
%Finding sigma
if problem==1 || problem==4 %no exterior term
    sigma = 0*ones(size(theta));
else
    sigma = sigsolve(theta,Pe,n,eps,randFun);
end

%Finding omega
if problem==0 || problem==3 %no interior term
    omega = 0*ones(size(theta));
else
    Fd = (AAA_LS_solve(z.',N,alpha,centr)).'; 
    omega = real(zedze.*(Fd-alpha*conj(z)/2));
end

%Finding C
if problem==0 || problem==1 || problem==2 %no area conservation
    C=0;
else 
    intsig = trapz(theta,sigma); intom = trapz(theta,omega); %integrals of sigma and omega
    intfz = integral(@(alpha) adzeta(alpha,c,N),-pi,pi); %integral of |f_zeta|
    C = -(intsig+beta*intom)/intfz;
end
end