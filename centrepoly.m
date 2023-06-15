function zc = centrepoly(z)
% Centre Polygon Function
%   zc = centrepoly(z)
%   centres the polygon z at the origin.
%
% INPUTS
%   z           = polygon data. List of polygon boundary points z=x+iy with
%                   coordinates (x,y).
%
% OUTPUTS   
%   zc          = recentred polygon data.
%                   
% END OF DOCUMENTATION
%
%Code
xpts = real(z); ypts = imag(z); pts = [xpts;ypts].';
[xc,yc]=centroid(polyshape(pts)); polycentre = xc+1i*yc; %centre of polygon
zc = z-polycentre; %zc = polygon centred at origin
end