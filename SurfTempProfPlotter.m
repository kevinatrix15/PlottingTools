%% This script parses and plots surface temperatre profiles along scan-aligned axis directions (x', y')
clear all, clc, close all
tot_steps = 1800;

for t=74:tot_steps+1
    filename1 = ['Result_T', num2str(t-1), '_SurfTempProfiles.csv'];
    M = csvread(filename1, 1);
    
    xmax = max(M(:,1));
    xmin = min(M(:,1)); 
    ymax = max(M(:,2 ));
    ymin = min(M(:,2));
    
    x_end = find(M(:,1) == xmax)           % Ending index of x' values

    xx = M(1:x_end,1), xprof = M(1:x_end,4);
    figure(2)
    subplot(1,2,1)
    plot(xx, xprof, 'bo'), hold on
%     axis([-4e-4 4e-4 0 4000])
    xlabel('x'' [m]'), ylabel('T [C]');
    title(['t = ', num2str(t-1)])
    
    
    yy = M(x_end+1:end,2), yprof = M(x_end+1:end,4);
%     figure(3)
    subplot(1,2,2)
    plot(yy, yprof, 'ro'), hold on
    axis([-4e-4 4e-4 0 4000])
    xlabel('y'' [m]'), ylabel('T [C]');
    
    
    plot_tool_dir = 'C:\Users\Kevin\OneDrive\3DSIM'
    curr_dir = pwd;
    cd(plot_tool_dir)
    %     x0 = [4000; -4.0e-006; 1.45e-005; 380; -0.004; 0.001];
%     x0 = [max(xprof); 0; 1e-6; 1e-6; 0; 1e-6; min(xprof)];
    x0 = [max(xprof); 0; 3e-6; 3733.32; -0.00020238; 1e-6; min(xprof)];
    options = optimset('Algorithm','levenberg-marquardt');
    [X, resnorm, residual, exitflag, output] = lsqcurvefit(@Gauss2Eq1D_funFit, x0, xx, xprof, [], [], options);
    figure(2)
    subplot(1,2,1)
    Gauss2Eq1DPlotter(X(1), X(2), X(3), X(4), X(5), X(6), X(7), -4e-4, 4e-4, '-g')
    cd(curr_dir)
%     pause(0.5)
%     hold off
    %     figure(2)
%     plot(M(x_end+1:end,2), M(x_end+1:end,4), '-x')


        filename2 = ['Result_T', num2str(t-1), '_Gauss2_1DParams.csv'];
        N = csvread(filename2, 1);
        a1x = N(:,1);        b1x = N(:,2);        c1x = N(:,3);
        a2x = N(:,4);        b2x = N(:,5);        c2x = N(:,6);
        dx = N(:,7);
        
        a1y = N(:,8);        b1y = N(:,9);        c1y = N(:,10);
        a2y = N(:,11);       b2y = N(:,12);       c2y = N(:,13);
        dy = N(:,14);
        
        figure(2)
        subplot(1,2,1)
        Gauss2Eq1DPlotter(a1x, b1x, c1x, a2x, b2x, c2x, dx, -4e-4, 4e-4,'--k');
        axis([-4e-4 4e-4 0 5000])
        hold off
        
%         figure(3)
        subplot(1,2,2)
        Gauss2Eq1DPlotter(a1y, b1y, c1y, a2y, b2y, c2y, dy, -4e-4, 4e-4, '--k');
        axis([-4e-4 4e-4 0 5000])
        pause(0.2)
        hold off
end


%% Fit gaussian after removing liquid nodes
clear all, clc
T_liq = 1690.15;

t = 3;

filename1 = ['Result_T', num2str(t-1), '_SurfTempProfiles.csv'];
M = csvread(filename1, 1);

xmax = max(M(:,1));
xmin = min(M(:,1));
ymax = max(M(:,2));
ymin = min(M(:,2));

x_end = find(M(:,1) == xmax)           % Ending index of x' values

T = M(:,4);
% T_x = zeros(xmax, 1);
% x_loc = zeros(xmax, 1);

c = 1;
% for i=1:x_end
for i=x_end+1:length(M)
    if T(i) < T_liq
        y_loc(c) = M(i,2);
        T_x(c) = T(i);
        c = c+1;
    end
end

%%
        % y_idx = [find(M(:,1) == xmax)+1:size(M,1)]
        % X =
        
        filename2 = ['Result_T', num2str(t-1), '_Gauss2_1DParams.csv'];
        N = csvread(filename2, 1);
        a1x = N(:,1);
        b1x = N(:,2);
        c1x = N(:,3);
        a2x = N(:,4);
        b2x = N(:,5);
        c2x = N(:,6);
        dx = N(:,7);
        
        a1y = N(:,8);
        b1y = N(:,9);
        c1y = N(:,10);
        a2y = N(:,11);
        b2y = N(:,12);
        c2y = N(:,13);
        dy = N(:,14);
        figure(3)
        Gauss2Eq1DPlotter(a1x, b1x, c1x, a2x, b2x, c2x, dx);
        figure(4)
        Gauss2Eq1DPlotter(a1y, b1y, c1y, a2y, b2y, c2y, dy);
        %     hold on
        
        figure(1)
        subplot(1,2,1)
        plot(M(1:x_end,1), M(1:x_end,4), 'bo')
        axis([-4e-4 4e-4 0 7e4])
        xlabel('x'' [m]'), ylabel('T [C]');
        hold off
        
        subplot(1,2,2)
        plot(M(x_end+1:end,2), M(x_end+1:end,4), 'ro')
        axis([-4e-4 4e-4 0 7e4])
        xlabel('y'' [m]'), ylabel('T [C]');
        pause(0.25)
        
%%
        % y_idx = [find(M(:,1) == xmax)+1:size(M,1)]
        % X =
        
        filename2 = ['Result_T', num2str(t-1), '_SkewNormTempCoefs.csv'];
        N = csvread(filename2, 1);
        xi_x = N(:,1);
        omeg_x = N(:,2);
        alph_x = N(:,3);
        xi_y = N(:,4);
        omeg_y = N(:,5);
        alph_y = N(:,6);
        figure(3)
        SkewNormPlotterFunc(xi_x, omeg_x, alph_x);
        %     hold on
        
        figure(1)
        subplot(1,2,1)
        plot(M(1:x_end,1), M(1:x_end,4), 'bo')
        axis([-4e-4 4e-4 0 7e4])
        xlabel('x'' [m]'), ylabel('T [C]');
        hold off
        
        subplot(1,2,2)
        plot(M(x_end+1:end,2), M(x_end+1:end,4), 'ro')
        axis([-4e-4 4e-4 0 7e4])
        xlabel('y'' [m]'), ylabel('T [C]');
        pause(0.25)
        
