function penguin_plots(M,spacing,fignumber)
% Penguin Plotting Function
%   penguin_plots(M,fignumber)
%   plot of free boundary shapes of cell array M in figure(fignumber).
%
% INPUTS
%   M           = cell array of free boundary data.
%                   Each entry (t-1) gives a list of complex z=x+iy values 
%                   for points with coordinates (x,y) on the free boundary 
%                   at time step t. (t=0: initial shape.)
%
%   spacing     = e.g. if spacing=n, then every nth step is plotted.
%
%   fignumber   = figure number. E.g. = 1 will display results in fig 1.
%
% OUTPUTS   
%   figure(fignumber)   = plot of free boundary shapes of cell array M. 
%
% END OF DOCUMENTATION
%
%Code
plotend = size(M,2); 

zinit = M{1}; zfin = M{end}; %initial and final free boundary shapes
minx = -0.5+min(real(zinit)); maxx=0.5+max(real(zfin)); yrange = (abs(maxx)+abs(minx)); %x and y ranges
axlist = [minx maxx -yrange/2 yrange/2]; %axis list

figure(fignumber)
for k=1:spacing:plotend
    z=M{k}; plot(real(z),imag(z),'LineWidth',2.5), hold on %free boundary plots
end
hold off, 
set(gca,'XColor', 'none','YColor','none')
axis(axlist), daspect([1 1 1]),

end