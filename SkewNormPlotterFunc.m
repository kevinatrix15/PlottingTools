%% Skew Normal function- Plotter

function [] = SkewNormPlotterFunc(xi1, omeg1, alph1)
fx = @(x)1/(omeg1*sqrt(2*pi)) * exp(-0.5*((x-xi1)/omeg1)^2) * (1+erf(alph1/sqrt(2) * (x-xi1)/omeg1))+353

%     subplot(1,2,1)
xmin = -10e-4, xmax = 10e-4;
% xmin = -10, xmax = 10;
fplot(fx, [xmin xmax])
% axis([xmin xmax 0 60])
end