function [RMSE, AreaError] = penguin_error(M)
% Penguin Error Function
%   [RMSE, AreaError] = penguin_error(M)
%   returns the root mean-squared error and area error for the 
%   free boundary shapes in cell M. 
%
% INPUTS
%   M           = cell array of free boundary data.
%                   Each entry (t-1) gives a list of complex z=x+iy values 
%                   for points with coordinates (x,y) on the free boundary 
%                   at time step t. (t=0: initial shape.)
% OUTPUTS   
%
%   RMSE        = root mean-squared error. Calculate the square of the
%                   difference between the point z on the free boundary 
%                   between time steps t and t-1. Sum these squares, divide
%                   by the number of points and square root this sum. The
%                   sum approaches a constant k for large time, so create a
%                   relative error with RMSE_exact = k.
%
%   AreaError   = area error. Finds the relative error of the area, where
%                   area_exact = initial area.
%
% NOTE
%   RMSE and AreaError of interest for the full penguin problem only i.e.
%   problem = 5 -- see penguin_RUN. Other values of problem variable may
%   give unexpected behaviour for RMSE and AreaError.
%
% END OF DOCUMENTATION
%
%Code
tsteps = size(M,2); %total number of time steps
rmse = 0*ones(1,tsteps); area = 0*ones(1,tsteps);

for k=1:tsteps
    if k~=1
        z2 = M{k}; n=size(z2,2); totalsum=0; %n=number of points on polygon
        for j=1:n
            totalsum=totalsum+(abs(z2(j)-z(j))^2); %sum of square of difference between jth point at current and previous time step
        end
        rmse(k)=sqrt(totalsum./n); %total sum divided by number of points
    end
    z=M{k}; area(k)=polyarea(real(z),imag(z)); %area of polygon
end 

rexact = rmse(end); %exact value of RMSE - the constant value it reaches for large time
aexact = area(1); %exact value of area - its initial value

RMSEmat = abs((rmse-rexact)./rexact); RMSE = RMSEmat(2:end-1); %RMSE ill defined at t=0 and t=end
Areamat = abs((area-aexact)./aexact); AreaError = Areamat(2:end); %AreaError ill defined at t=0
end