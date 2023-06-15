function res = odeqns(t,c,cdot,N,n,Pe,beta,alpha,eps,randFun,problem)
% OD Equations Function
%   res = odeqns(t,c,cdot,N,n,Pe,beta,alpha,eps,randFun,problem)
%   represents the system of ODEs to be solved by decic/ode15i.
%   Recall the conformal map z = f(zeta,t).
%
% INPUTS
%   t           = time variable. (Unused in odeqns but needed for ode15i).
%   c           = initial Laurent coefficient.
%   cdot        = initial Laurent coefficient time derivative.
%   N           = series truncation of conformal map. 
%   n           = number of unknowns = 2N+3.
%   Pe          = Peclet number - strength of the flow.
%   beta        = thermal conductivity ratio - strength of interior effect.
%   alpha       = Poisson equation forcing.
%   eps         = amplitude of noise.
%   randFun     = random (spatial) function: assigns each position theta a 
%                   value between -1 and 1.
%   problem     = problem choice. The ode solver is applicable to six free
%                   boundary problems depending on the inclusion/exclusion
%                   of the exterior, interior and area conserving terms in
%                   the Stefan boundary condition. See [1] for more
%                   details.
%                    0 (ext); 1 (int); 2 (ext+int); 3 (ext+areaC), 
%                       4 (int+areaC), 5 (ext+int+areaC).
%
% OUTPUTS
%   res         = output used by decic/ode15i.
%
% REFERENCES
%   [1]	        Harris, S.J., McDonald, N.R. (2023) "Penguin Huddling: 
%               A Continuum Model". Acta Appl. Math. 185, 7. 
%               https://doi.org/10.1007/s10440-023-00578-2.
%
% END OF DOCUMENTATION
%
%Code
theta = linspace(0,2*pi,n+1); zeta=exp(1i*theta); %end point theta=2*pi is not included
centr = c(2)+1i*c(N+3); %conformal centre

z=c(1).*zeta; %z
dzeta = c(1); %dz/dzeta
dt=cdot(1).*zeta; %dz/dt
for m=2:N+2
    z=z+(c(m)+1i*c(N+1+m))*zeta.^(-(m-2));
    dzeta=dzeta-(m-2).*(c(m)+1i*c(N+1+m))*zeta.^(-(m-2)-1);
    dt=dt+(cdot(m)+1i*cdot(N+1+m))*zeta.^(-(m-2));
end
zedze = zeta.*dzeta; %zeta.*dz/dzeta
[sigma,omega,C] = SOCcalc(theta,z,zedze,Pe,beta,alpha,eps,randFun,c,N,n,centr,problem); %sigma, omega, C calculations

for k=1:n
res(k)=(-Pe.*real(dt(k)*conj(zedze(k))))-sigma(k)-beta.*omega(k)-C*abs(zedze(k)); res=res.'; %Polubarinova-Galin equation
end
res=res.'; %output used in decic/ode15i
end