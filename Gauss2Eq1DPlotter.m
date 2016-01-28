%% 2 Eqn Gaussian function- Plotter

function [] = Gauss2Eq1DPlotter(a1, b1, c1, a2, b2, c2, d, xmin, xmax, lstyle)
fx = @(x)a1*exp(-((x-b1)/c1)^2) + a2*exp(-((x-b2)/c2)^2) + d;

f1 = @(x)a1*exp(-((x-b1)/c1)^2) + d;
f2 = @(x)a2*exp(-((x-b2)/c2)^2) + d;
%     subplot(1,2,1)
% xmin = -4e-6, xmax = 4e-6;
% xmin = -10, xmax = 10;
% figure(2)
% subplot(1,2,1)
fplot(fx, [xmin xmax], lstyle)
hold on
fplot(f1, [xmin xmax], '--g')
fplot(f2, [xmin xmax], '--m')
plot([0 b1],[a1+d a1+d], '-r',[b1 b1],[0+d a1+d], '-r',[0 b2],[a2+d a2+d], '-b',[b2 b2],[0+d a2+d], '-b')

% h = ezplot(fx, [xmin xmax]);
% set(h, 'LineSpec',lstyle,'LineWidth',2);

% axis([-4e-4 4e-4 0 4000])
% axis([xmin xmax 0 60])
end