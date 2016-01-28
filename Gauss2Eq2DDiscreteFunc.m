function [hxy] = Gauss2Eq2DDiscreteFunc(p, P0, th, x, y)

xp = (x - P0(1))*cos(th) + (y-P0(2))*sin(th);
yp = (P0(1) - x)*sin(th) + (y-P0(2))*cos(th);

% xp = x; yp = y;

hxy = (p(1)*exp(-((xp-p(2))/p(3))^2 - ((yp-p(8))/p(9))^2) ...
    + p(10)*exp(-((xp-p(2))/p(3))^2 - ((yp-p(11))/p(12))^2) ...
    + p(7)*exp(-((xp-p(5))/p(6))^2 - ((yp-p(8))/p(9))^2) ...
    + p(4)*exp(-((xp-p(5))/p(6))^2 - ((yp-p(11))/p(12))^2) ...
    + p(13));


end