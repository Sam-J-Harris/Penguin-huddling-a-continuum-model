function penguin_plots_heat_flux(M,Cflux,Pflux,fignumber1,fignumber2,fignumber3)
% Penguin Heat Flux Plotting Function
%   penguin_plots_heat_flux(M,Cflux,Pflux,fignumber1,fignumber2,fignumber3)
%   three plots of heat flux data across the steady shape of free boundary 
%   data in cell array M.
%
% INPUTS
%   M           = cell array of free boundary data.
%                   Each entry (t-1) gives a list of complex z=x+iy values 
%                   for points with coordinates (x,y) on the free boundary 
%                   at time step t. (t=0: initial shape.)
%
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
%   fignumber1  = figure number. E.g. = 1 will display results in fig 1.
%
%   fignumber2  = ""
%
%   fignumber3  = ""
%
% OUTPUTS   
%   figure(fignumber1)   = plot of the upper half of the steady shape of
%                           cell array M. Also plotted are arrows normal to
%                           points on the boundary scaled by the physical
%                           heat flux.
%
%   figure(fignumber2)   = plot of canonical heat flux vs theta.
%
%   figure(fignumber3)   = plot of physical heat flux vs x.
%                   
% END OF DOCUMENTATION
%
%Code
z=centrepoly(M{end}); x1=real(z); y1=imag(z); %polygon centred at origin
theta=Cflux{1}; cflux=Cflux{2}; z2=Pflux{1}; pflux=Pflux{2};

dztemp = z2 - circshift(z2,-1); v = 1i.*(dztemp)./abs(dztemp); %used in calculating normal vector 
nvec = abs(pflux).*(v+circshift(v,1))./abs(v+circshift(v,1)); %outward unit normal vector scaled by pflux

x2=real(z2); y2=imag(z2); u = real(nvec); v = imag(nvec); %inputs for quiver function

figure(fignumber1)
plot(x1,y1,'LineWidth',2.5), hold on, quiver(x2(1:2:end),y2(1:2:end),u(1:2:end),v(1:2:end),'LineWidth',2), hold off,
set(gca,'YColor','none')
daspect([1 1 1]), axis([-1.75 1.75 0 1.75]), axis square, daspect([1 1 1]),

figure(fignumber2)
plot(theta(1:end/2),cflux(1:end/2),'LineWidth',2.5),daspect([1 1 1]), axis square,
set(gca,'XTick',0:pi/4:pi) 
set(gca,'XTickLabel',{'0','\pi/4','\pi/2','3\pi/4','\pi'},'FontSize',15)

figure(fignumber3)
plot(x2(1:end/2),pflux(1:end/2),'LineWidth',2.5), daspect([1 1 1]), axis square,

end