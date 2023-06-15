function penguin_plots_steady_shape(totalM,fignumber)
% Penguin Steady Shape Plotting Function
%   penguin_plots_steady_shape(totalM,fignumber)
%   plots of the steady shapes of free boundary data in cell array totalM 
%   in figure(fignumber). Steady shapes are recentred to the origin. For
%   definition of steady shape, see [1].
%
% INPUTS
%   totalM      = cell array of cell arrays of free boundary data.
%                   E.g. totalM = {M1, M2, M3}, where each cell array e.g.
%                   M1 is a cell array of free boundary data.
%
%   fignumber   = figure number. E.g. = 1 will display results in fig 1.
%
% OUTPUTS   
%   figure(fignumber)   = plot of steady shapes recentred at the origin 
%                           of free boundary data in cell array totalM. 
%
% REFERENCES
%   [1]	        Harris, S.J., McDonald, N.R. (2023) "Penguin Huddling: 
%               A Continuum Model". Acta Appl. Math. 185, 7. 
%               https://doi.org/10.1007/s10440-023-00578-2.
%
% END OF DOCUMENTATION
%
%Code
plotend = size(totalM,2); 

figure(fignumber)
for k=1:plotend
    z=totalM{k}; zc = centrepoly(z{end}); %zc = recentred polygon
    plot(real(zc),imag(zc),'LineWidth',2), hold on %steady shape plots
end
hold off, 
set(gca,'XColor', 'none','YColor','none')
daspect([1 1 1]), axis([-1.5 1.5 -1.5 1.5]), 
if plotend>1
    legend()
end

end