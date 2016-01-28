function ellipse(x0, y0, a, b, lcolor)

% a = horizontal radius, b = vertical radius, x0,y0 = coordinates for
% center of ellipse

t=-pi:0.01:pi;
x=x0+a*cos(t);
y=y0+b*sin(t);
plot(x,y, lcolor)