%% Fits Gaus2 function to experimental temperature data
function fx = Gauss2Eq1D_funFit(x, xdata)

fx = x(1)*exp(-((xdata-x(2))/x(3)).^2) + x(4)*exp(-((xdata-x(5))/x(6)).^2) + x(7);
