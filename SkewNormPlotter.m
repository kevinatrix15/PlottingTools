%% Skew Normal function- Plotter

% xi1 = -0.0017173113660474463;  omeg1 = 0.027199098075252996;  alph1 = 23.464263673642549
% xi1 = 0.0056;  omeg1 = 0.073;  alph1 = 0.59
% a = 4014;   b = -0.000115;  c = 0.0001445;  alph1 = -0.0004366;     d = 353;
a = 3800;   b = -0.000115;  c = 0.0001445;  alph1 = 0.04366;     d = 353;

% 
% a2 = 0.15;  b2 = -0.25;
% a3 = 0.15;  b3 = -0.125;
% a4 = 0.15;  b4 = -0.15;
% mu = 0.2;

% fx = @(x)1/(omeg1*sqrt(2*pi)) * exp(-0.5*((x-xi1)/omeg1)^2) * (1+erf(alph1/sqrt(2) * (x-xi1)/omeg1))
fx = @(x)a * exp(-0.5*((x-b)/c)^2) * (1+erf(alph1/sqrt(2) * (x-b)/c)) + d

% gx = @(x)a2*exp(-(x-mu).^2/(2*b2^2))
% fy = @(y)a3*exp(-(y-mu).^2/(2*b3^2))
% gy = @(y)a4*exp(-(y-mu).^2/(2*b4^2))


% mxy = @(x,y) hx(x)*hy(y)
% ezplot(hx, hy, [-9, 9]), hold on
% ezplot(mxy), hold on
figure(5)
fplot(fx, [-5e-4 5e-4]), hold on
% fplot(fy, [0 1])
% fplot(gy, [0 1])