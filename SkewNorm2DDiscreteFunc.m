function [hxy] = SkewNorm2DDiscreteFunc(p, P0, th, x, y)

xp = (x - P0(1))*cos(th) + (y-P0(2))*sin(th);
yp = (P0(1) - x)*sin(th) + (y-P0(2))*cos(th);

ax = p(1);  bx = p(2);  cx = p(3);  alphx = p(4);
ay = p(5);  by = p(6);  cy = p(7);  alphy = p(8);
d = p(9);

% if (abs(ax-ay) < 10)
    a = 0.5*(ax+ay);
    
erfx = erf(alphx/sqrt(2)*(xp-bx)/cx);
erfy = erf(alphy/sqrt(2)*(yp-by)/cy);
hxy = a*exp(-0.5*((xp-bx)/cx)^2 - 0.5*((yp-by)/cy)^2)*...
    (1 + erfx + erfy + erfx*erfy) + d;

end