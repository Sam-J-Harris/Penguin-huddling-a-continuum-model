%% Penguin Free Boundary Problem RUN
%   RUN this code to obtain the evolution of the free boundary.
%   Change the below inputs as desired.
%  
% END OF DOCUMENTATION
%
%Code
clc; close all; clear; set(0,'DefaultFigureVisible','on'); warning('on','all');

%Inputs
Pe = 10.0;          %Peclet number.
beta = 0.1;         %thermal conductivity ratio.
alpha = 3;          %Poisson eq. forcing.
eps = 0;            %amplitude of noise -- set to 0 for no noise
N=128;              %series truncation of conformal map. Number of unknowns = n = 2N+3
steps = 20;         %number of time steps.
tstep = 0.1*Pe;     %value of each time step.
toly = 1e-8;        %relative and absolute tolerances on ode15i.
problem=5;          %problem choice: exterior, interior, area conserving terms. 0 (ext); 1 (int); 2 (ext+int); 3 (ext+areaC), 4 (int+areaC), 5 (ext+int+areaC).
shape = 21;         %initial shape: 21 (circle), 22 (slanted ellipse), 24 (triangle), 51 (irregular pentagon) -- see penguin_initial_shape for more.

% Free boundary evolution solver and plots
[M,Cflux,Pflux] = penguin_ode_solve(Pe,beta,alpha,eps,N,steps,tstep,toly,problem,shape); % M is a 1x(steps+1) cell of the free boundary shapes; cell entry i is the free boundary at step (i-1) -- step 0 is the initial shape.
[RMSE, AreaError] = penguin_error(M); %finds the root mean-squared error and area error -- see penguin_error for more detail.

penguin_plots(M,4,1); %plot of free boundary propagation at every time step (spacing=1) in figure 1.
penguin_plots_steady_shape({M},2); %plots the steady shape in figure 2. First input must be a cell array of cell arrays.
penguin_plots_heat_flux(M,Cflux,Pflux,3,4,5); %heat flux plots in figures 3, 4 and 5.