function [] = SkewNorm1DPlotter(a, b, c, d, alph, xmin, xmax, lstyle)
fx = @(x) a * exp(-0.5*((x-b)/c)^2) * (1+erf(alph/sqrt(2) * (x-b)/c))+d

fplot(fx, [xmin xmax], lstyle)
hold on

% plot([0 b1],[a1+d a1+d], '-r',[b1 b1],[0+d a1+d], '-r',[0 b2],[a2+d a2+d], '-b',[b2 b2],[0+d a2+d], '-b')
