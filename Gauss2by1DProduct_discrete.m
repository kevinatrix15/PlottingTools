function [hxy] = Gauss2by1DProduct_discrete(px, py, P0, th, x, y)

xp = (x - P0(1))*cos(th) + (y-P0(2))*sin(th);
yp = (P0(1) - x)*sin(th) + (y-P0(2))*cos(th);


fx = px(1)*exp(-((xp-px(2))/px(3))^2) + px(4)*exp(-((xp-px(5))/px(6))^2);
fy = py(1)*exp(-((yp-py(2))/py(3))^2) + py(4)*exp(-((yp-py(5))/py(6))^2);

d = 0.5*(px(7)+py(7));

hxy = fx * fy + d;

end