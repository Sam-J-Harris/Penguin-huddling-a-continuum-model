function absdzeta = adzeta(alpha,c,N)
% Finding |dz/dzeta| Function
%   absdzeta = adzeta(alpha,c,N)
%   function defined in order to use "integrate" function -- see SOCcalc.
%   (Note: alpha used to avoid confusion with theta in SOCcalc)
%
% INPUTS
%   alpha       = azimuthal angle of points around the zeta-disk.
%   c           = initial Laurent coefficient.
%   N           = series truncation of conformal map. 
%
% OUTPUTS
%   absdzeta    = |dz/dzeta|. Recall that z = f(zeta,t).
%
% END OF DOCUMENTATION
%
%Code
zeta = exp(1i*alpha); 
dzeta = c(1); %dz/dzeta
for m =2:N+2
    dzeta=dzeta+(-(m-2))*(c(m)+1i*c(N+1+m))*zeta.^(-(m-1));
end
absdzeta=abs(dzeta);
end