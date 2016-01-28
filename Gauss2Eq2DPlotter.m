%% 2 Eqn Gaussian function- Plotter

function [] = Gauss2Eq2DPlotter(p, P0, th, xmin, xmax, ymin, ymax)
% fx = @(x)p(1)*exp(-((x-p(5))/p(6))^2) + p(2)*exp(-((x-p(9))/p(10))^2);
% fy = @(y)p(3)*exp(-((y-p(7))/p(7))^2) + p(4)*exp(-((y-p(11))/p(12))^2);
% hxy = @(x,y)(fx(x).*fy(y) + p(13));
syms x y xp yp
th = -th;
% Rotate and translate coordinates for b-params
% p(2) = p(2)*cos(-th) + P0(1);
% p(5) = p(5)*cos(-th) + P0(1);
% p(8) = p(8)*cos(-th) + P0(2);
% p(11) = p(11)*cos(-th) + P0(2);
p(2) = (p(2) - P0(1))*cos(th) - P0(2)*sin(th);
p(5) = (p(5) - P0(1))*cos(th) - P0(2)*sin(th);
p(8) = P0(1)*sin(th) + (p(8) - P0(2))*cos(th);
p(11) = P0(1)*sin(th) + (p(11) - P0(2))*cos(th);

hxy = @(xp,yp)(p(1)*exp(-((xp-p(2))/p(3)).^2 - ((yp-p(8))/p(9)).^2) ...
    + p(10)*exp(-((xp-p(2))/p(3)).^2 - ((yp-p(11))/p(12)).^2) ...
    + p(7)*exp(-((xp-p(5))/p(6)).^2 - ((yp-p(8))/p(9)).^2) ...
    + p(4)*exp(-((xp-p(5))/p(6)).^2 - ((yp-p(11))/p(12)).^2) ...
    + p(13));

xp_n = (x - P0(1))*cos(th) + (y-P0(2))*sin(th);
yp_n = (P0(1) - x)*sin(th) + (y-P0(2))*cos(th);
subs(xp, xp_n);
subs(yp, yp_n);
% hxy = @(x,y)(p(1)*exp(-((x-p(5))/p(6)).^2 - ((y-p(9))/p(10)).^2) ...
%             + p(2)*exp(-((x-p(5))/p(6)).^2 - ((y-p(11))/p(12)).^2) ...
%             + p(3)*exp(-((x-p(7))/p(8)).^2 - ((y-p(9))/p(10)).^2) ...
%             + p(4)*exp(-((x-p(7))/p(8)).^2 - ((y-p(11))/p(12)).^2) ...
%             + p(13));
% hxy = @(x,y)(p(1)*exp(-((x-p(3))/p(4)).^2 - ((y-p(7))/p(8)).^2) ...
%             + p(2)*exp(-((x-p(5))/p(6)).^2 - ((y-p(7))/p(8)).^2) ...
%             + p(9));
% colormap(jet);
% contourf(qx,qy,qz, [350 400 500 600 700 1000 1200 1500 1800 1923 3315 5000 8000 10000]), hold on;
% caxis([350, 3500])
% caxis manual;
% plot(SallyT(:,1), SallyT(:,2), '*'), hold off

% Map analytical solution to discrete grid for plotting 2D contourf:
n = 150;
[X,Y] = meshgrid(0:1:n-1, 0:1:n-1);
X = 

R = [xmin, xmax, ymin, ymax];
ezsurf(hxy, R, 90)


end