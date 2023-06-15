function Fd=AAA_LS_solve(Z,N,alpha,centr)
% AAA-LS Function
%   Fd=AAA_LS_solve(Z,N,alpha,centr)
%   returns F'(z) where u=Re(F) is a harmonic function in the interior -- see 
%   [1] for conditions on u. Uses AAA and least squares algorithms -- see
%   e.g [2] and [3]. Chebfun package must be installed -- see [4].
%
% INPUTS
%   Z           = list of complex z=x+iy values representing points with 
%                   coordinates (x,y) on the free boundary.
%   N           = series truncation of conformal map. 
%   alpha       = Poisson equation forcing.
%   centr       = (conformal) centre of the shape traced by the free boundary.
%
% OUTPUTS
%   Fd          = F'(z) where u=Re(F) is harmonic in the domain.
%
% REFERENCES
%   [1]	        Harris, S.J., McDonald, N.R. (2023) "Penguin Huddling: 
%               A Continuum Model". Acta Appl. Math. 185, 7. 
%               https://doi.org/10.1007/s10440-023-00578-2.
%
%   [2]         Costa, S., & Trefethen, L. N. (2021). "AAA-least squares 
%               rational approximation and solution of Laplace problems". 
%               arXiv preprint arXiv:2107.01574.
%
%   [3]         Nakatsukasa Y., Sete 0., Trefethen L.N. (2018) "The AAA algorithm
%               for rational approximation". SIAM J. Sci. Comp. 40, A1494-A1522.
%
%   [4]         Driscoll, T. A., Hale, N., & Trefethen, L. N. (2014). "Chebfun guide".
%               Pafnuty Press, Oxford, 2014; see also www.chebfun.org
%
% END OF DOCUMENTATION
%
%Code
Na=N+2; Za = circshift(Z,N+1); %points around polygon shifted to avoid ill behaviour at the start/end points
h = @(z) 1+alpha*0.25.*z.*conj(z); H = h(Za); %boundary condition
inpolygonc = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w)); %determines if point is in polygon

% Global AAA poles %%%%%%%
[~,polk] = aaa(H,Za,'cleanup',1,'toly',1e-8); pol = (polk(~inpolygonc(polk,Za))).'; %pol = poles outside polygon
if size(pol,1)==0 && size(pol,2)==0 %if there are no poles, set pol=d=0
    pol = 0; d = 0;
else
    d = min(abs((Za)-pol),[],1); %d=distance between pole and closest point on polygon
end

% Least-Squares %%%%%%%% 
[Hes,P] = VAorthog(Za-centr,Na); Q = d./((Za)-pol); %Q scaled by d
A = [real(P) real(Q) -imag(P) -imag(Q)]; c = reshape(A\H,[],2)*[1; 1i];  
fd = @(z) reshape([VAevald((z(:)-centr),Hes) -d./(((z(:)-pol)).^2)]*c,size(z)); %note the function VAevald, gives derivative output only from VAeval
Fd = fd(Z); %gives F'(z) where u = Re(F).
end