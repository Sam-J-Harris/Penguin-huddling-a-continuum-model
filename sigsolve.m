function sigma = sigsolve(theta,Pe,n,eps,randFun) 
% Sigma Solver Function
%   sigma = sigsolve(theta,Pe,n,eps,randFun) 
%   solves the integral equation governing sigma, see Appendix A of [1].
%
% INPUTS
%   theta       = azimuthal angle of points around the zeta-disk.
%   Pe          = Peclet number - strength of the flow.
%   n           = number of unknowns (in conformal map f).
%   eps         = amplitude of noise.
%   randFun     = random (spatial) function: assigns each position theta a 
%                   value between -1 and 1.
%
% OUTPUTS
%   sigma       = exterior effect term in the PG equation.
%
% REFERENCES
%   [1]         Ladd, A. J., Yu, L., & Szymczak, P. (2020). "Dissolution of 
%               a cylindrical disk in Hele-Shaw flow: a conformal-mapping 
%               approach". J. Fluid Mech., 903.
%
% END OF DOCUMENTATION
%
%Code
Nt = (n-1)/2; Kmat = 0*ones(Nt,Nt); %sigma in range (0,pi), endpoints not included
deltaK = abs((theta(2)-theta(1))); %deltaK = spacing between theta values
for m = 1:Nt
    K = @(ptheta) exp(Pe.*(cos(theta(m+1))-cos(ptheta))).*besselk(0, abs(Pe.*(cos(theta(m+1))-cos(ptheta)))); %kernel function -- see [1] (A2), (A4)
    Kmat(m,:) = deltaK*K(theta(2:Nt+1).'); Kmat(m,m) = 0; %non-diagonal elements -- see [1] (A4)
    Ksum = sum(Kmat,2);  Kmn = Ksum(m); 
    Im = integral(K,0,theta(m+1))+integral(K,theta(m+1),pi); %integral of K between 0 and pi, care taken at theta(m+1) point
    Kmat(m,m) = Im - Kmn; %diagonal elements -- see [1] (A6)
end
pivec = (pi*ones(1,Nt)).'; sigmain = (Kmat\pivec).'; %main sigma calculation -- see [1] (A5)
sigmacal = [1/pi sigmain flip(sigmain) 1/pi]; %sigma values for [0,2*pi], sigma(0)=sigma(2*pi)=1/pi and sigma(pi+theta)=sigma(theta)
sigma = (1+eps.*randFun).*sigmacal; %added noise to sigma (optional) -- eps=0 is the standard case
end