function penguin_plots_aspect_ratio(totalM1,totalM2,totalM3,parlist,fignumber)
% Penguin Aspect Ratio Plotting Function
%   penguin_plots_aspect_ratio(totalM1,totalM2,totalM3,parlist,fignumber)
%   plots the aspect ratio (AR) of the final steady shape of each
%   experiment in totalM1, totalM2 and totalM3 vs the parameter list
%   specified in parlist.
%
% INPUTS
%   totalM1      = cell array of cell arrays of free boundary data.
%                   E.g. totalM = {M1, M2, M3}, where each cell array e.g.
%                   M1 is a cell array of free boundary data.
%
%   totalM2      = ""
%
%   totalM3      = ""
%
%   parlist      = list of parameters (either Pe or beta).
%
%   fignumber   = figure number. E.g. = 1 will display results in fig 1.
%
% OUTPUTS   
%   figure(fignumber)   = plot of AR vs parlist: three curves for the three
%                           data inputs totalM1, totalM2, totalM3.
%
% REFERENCES
%   [1]	        Harris, S.J., McDonald, N.R. (2023) "Penguin Huddling: 
%               A Continuum Model". Acta Appl. Math. 185, 7. 
%               https://doi.org/10.1007/s10440-023-00578-2.
%
% END OF DOCUMENTATION
%
%Code
plotend = size(totalM1,2); arlist1=[]; arlist2=[]; arlist3=[];

for k=1:plotend %compute AR from lists of free boundary shapes
    z1 = totalM1{k}; zc1 = centrepoly(z1{end}); sz1 =size(zc1,2); 
    x1 = real(zc1); y1 = imag(zc1); % centered steady shape of totalM1{k}
    a1 = abs(x1(1)-x1(sz1/2)); b1 = abs(y1(sz1/4)-y(3*sz1/4)); 
    ar1 = a1/b1; arlist1 = [arlist1; ar1]; % AR of totalM1{k} appended to list

    z2 = totalM2{k}; zc2 = centrepoly(z2{end}); sz2 =size(zc2,2); 
    x2 = real(zc2); y2 = imag(zc2); % centered steady shape of totalM2{k}
    a2 = abs(x2(1)-x2(sz2/2)); b2 = abs(y2(sz2/4)-y(3*sz2/4)); 
    ar2 = a2/b2; arlist2 = [arlist2; ar2]; % AR of totalM2{k} appended to list

    z3 = totalM3{k}; zc3 = centrepoly(z3{end}); sz3 =size(zc3,2); 
    x3 = real(zc3); y3 = imag(zc3); % centered steady shape of totalM3{k}
    a3 = abs(x3(1)-x3(sz3/2)); b3 = abs(y3(sz3/4)-y3(3*sz3/4)); 
    ar3 = a3/b3; arlist3 = [arlist3; ar3]; % AR of totalM3{k} appended to list
end

X=parlist; Y1=arlist1; Y2=arlist2; Y3=arlist3;

figure(fignumber)
plot(X,Y1,'-o','Color','k'), hold on,
plot(X,Y2,'-x','Color','b'),
plot(X,Y3,'-*','Color','r'),
hold off, 
if parlist(1)==0    % AR vs beta
    xlabel('$\beta$','FontSize',19,'interpreter','latex'),  lgd = {'Pe $=10$', 'Pe $=50$', 'Pe $=100$'}; lc = 'northeast';
else                % AR vs Pe
    xlabel('Pe','FontSize',19,'interpreter','latex'), lgd = {'$\beta=0$', '$\beta=0.1$', '$\beta=1$'}; lc = 'southeast';
end
ylabel('Aspect ratio','FontSize',19), legend(lgd,'interpreter','latex','location',lc,'FontSize',19), axis square
end