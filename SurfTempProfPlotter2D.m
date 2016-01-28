%% This script parses and plots surface temperatre profiles along scan-aligned axis directions (x', y')
clear all, clc, close all
tot_steps = 249;

for t=250:tot_steps+1
%     filename1 = ['Result_T', num2str(t-1), '_SurfTempProfiles.csv'];
%     M = csvread(filename1, 1);
%     
%     xmax = max(M(:,1));
%     xmin = min(M(:,1));
%     ymax = max(M(:,2 ));
%     ymin = min(M(:,2));
%     
%     x_end = find(M(:,1) == xmax)           % Ending index of x' values
%     
%     xx = M(1:x_end,1), xprof = M(1:x_end,4);
%     figure(2)
%     subplot(1,2,1)
%     plot(xx, xprof, 'bo'), hold on
%     %     axis([-4e-4 4e-4 0 4000])
%     xlabel('x'' [m]'), ylabel('T [C]');
%     title(['t = ', num2str(t-1)])
%     
%     
%     yy = M(x_end+1:end,2), yprof = M(x_end+1:end,4);
%     %     figure(3)
%     subplot(1,2,2)
%     plot(yy, yprof, 'ro'), hold on
%     axis([-4e-4 4e-4 0 4000])
%     xlabel('y'' [m]'), ylabel('T [C]');
    
    
    plot_tool_dir = 'C:\Users\Kevin\OneDrive\3DSIM'
%     curr_dir = pwd;
%     cd(plot_tool_dir)
%     %     x0 = [4000; -4.0e-006; 1.45e-005; 380; -0.004; 0.001];
%     %     x0 = [max(xprof); 0; 1e-6; 1e-6; 0; 1e-6; min(xprof)];
%     x0 = [max(xprof); 0; 3e-6; 3733.32; -0.00020238; 1e-6; min(xprof)];
%     options = optimset('Algorithm','levenberg-marquardt');
%     [X, resnorm, residual, exitflag, output] = lsqcurvefit(@Gauss2Eq1D_funFit, x0, xx, xprof, [], [], options);
%     figure(2)
%     subplot(1,2,1)
%     Gauss2Eq1DPlotter(X(1), X(2), X(3), X(4), X(5), X(6), X(7), -4e-4, 4e-4, '-g')
%     cd(curr_dir)
%     %     pause(0.5)
%     %     hold off
%     %     figure(2)
%     %     plot(M(x_end+1:end,2), M(x_end+1:end,4), '-x')
    
    
    filename2 = ['Result_T', num2str(t-1), '_Gauss_2DParams.csv'];
    N = csvread(filename2, 1);
    a1x = N(:,1);        b1x = N(:,2);        c1x = N(:,3);
    a2x = N(:,4);        b2x = N(:,5);        c2x = N(:,6);
   
    a1y = N(:,7);        b1y = N(:,8);        c1y = N(:,9);
    a2y = N(:,10);       b2y = N(:,11);       c2y = N(:,12);
    d = N(:,13);
    P0(1) = N(:,14);    P0(2) = N(:,15);    theta = N(:,16);
    
    p = [a1x b1x c1x a2x b2x c2x a1y b1y c1y a2y b2y c2y d];
    
%     xmin = -3e-4;   xmax = 3e-4;
    xmin = 0;   xmax = 6e-4;

    ymin = xmin;    ymax = xmax;
    curr_dir = pwd;
    cd(plot_tool_dir)
    figure(2)
%     subplot(1,2,1)
    Gauss2Eq2DPlotter(p, P0, theta, xmin, xmax, ymin, ymax);
    cd(curr_dir)
%     axis([0 1e-3 0 1e-3 -200 8000])
    axis([xmin xmax ymin ymax -200 8000])
    colorbar
%     hold off
%     
%     %         figure(3)
%     subplot(1,2,2)
%     Gauss2Eq1DPlotter(a1y, b1y, c1y, a2y, b2y, c2y, dy, -4e-4, 4e-4, '--k');
%     axis([-4e-4 4e-4 0 5000])
    pause
%     hold off
end

%% This script parses and plots surface temperatre profiles along scan-aligned axis directions (x', y')
clear all, clc, close all
tot_steps = 1800;
curr_dir = pwd;

%----- The below is for SallyReader -----
Nodes=26;
sd=8;
SallyT=csvread('Result_L1_O1_T0_Temperature.csv');
% set(gca, 'nextplot', 'replacechildren');
bedwidth_x=max(SallyT(:,1));
bedwidth_y=max(SallyT(:,2));
% figure
[qx,qy] = meshgrid([0:bedwidth_x/(sd*(Nodes-1)):bedwidth_x],[0:bedwidth_y/(sd*(Nodes-1)):bedwidth_y]);
csvfiles = dir('*Temperature.csv'); 
[~,index] = sortrows([csvfiles.datenum].'); csvfiles = csvfiles(index);
numfiles = length(csvfiles);

for t=1:1:tot_steps+1
    
    %----- 1D 2 term Ploting -----
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
    axis([-6e-4 6e-4 0 5000])
%     axis([-4e-4 4e-4 0 4000])
    xlabel('x'' [m]'), ylabel('T [K]');
    title(['t = ', num2str(t-1)])
    
    
    yy = M(x_end+1:end,2), yprof = M(x_end+1:end,4);
    subplot(1,2,2)
    plot(yy, yprof, 'ro'), hold on
    axis([-6e-4 6e-4 0 5000])
    %axis([-4e-4 4e-4 0 4000])
    xlabel('y'' [m]'), ylabel('T [K]');
    
    
    plot_tool_dir = 'C:\Users\Kevin\OneDrive\3DSIM'
    cd(plot_tool_dir)
    %     x0 = [4000; -4.0e-006; 1.45e-005; 380; -0.004; 0.001];
%     x0 = [max(xprof); 0; 1e-6; 1e-6; 0; 1e-6; min(xprof)];
    x0 = [max(xprof); 0; 3e-6; 3733.32; -0.00020238; 1e-6; min(xprof)];
    options = optimset('Algorithm','levenberg-marquardt');
    [X, resnorm, residual, exitflag, output] = lsqcurvefit(@Gauss2Eq1D_funFit, x0, xx, xprof, [], [], options);
    figure(2)
    subplot(1,2,1)
%     Gauss2Eq1DPlotter(X(1), X(2), X(3), X(4), X(5), X(6), X(7), -4e-4, 4e-4, '-g')
    cd(curr_dir)
    axis([-6e-4 6e-4 0 5000])
    %     pause(0.5)
    %     hold off
    %     figure(2)
    %     plot(M(x_end+1:end,2), M(x_end+1:end,4), '-x')
    
    
    filename2 = ['Result_T', num2str(t-1), '_Gauss2_1DParams.csv'];
    N = csvread(filename2, 1);
    a1x1d = N(:,1);        b1x1d = N(:,2);        c1x1d = N(:,3);
    a2x1d = N(:,4);        b2x1d = N(:,5);        c2x1d = N(:,6);
    dxx1d = N(:,7);
    
    a1y1d = N(:,8);        b1y1d = N(:,9);        c1y1d = N(:,10);
    a2y1d = N(:,11);       b2y1d = N(:,12);       c2y1d = N(:,13);
    dyy1d = N(:,14);
    
    P0(1) = N(:,15);    P0(2) = N(:,16);    theta = N(:,17);
    
    figure(2)
    subplot(1,2,1)
    Gauss2Eq1DPlotter(a1x1d, b1x1d, c1x1d, a2x1d, b2x1d, c2x1d, dxx1d, -6e-4, 6e-4,'--k');
    axis([-4e-4 2e-4 0 5000])
    legend('T_s_i_m', 'T_f_i_t')
    hold off
    
    subplot(1,2,2)
    Gauss2Eq1DPlotter(a1y1d, b1y1d, c1y1d, a2y1d, b2y1d, c2y1d, dyy1d, -6e-4, 6e-4, '--k');
    axis([-3e-4 3e-4 0 5000])
    %pause(0.2)
    hold off
    
    
    
%     %----- Product of 2 1D functions -----
%     filename3 = ['Result_T', num2str(t-1), '_Gauss_2DParams.csv'];
%     NN = csvread(filename3, 1);
%     P0(1) = NN(:,14);    P0(2) = NN(:,15);    theta = NN(:,16);
%     
%     px = [a1x b1x c1x a2x b2x c2x dxx];
%     py = [a1y b1y c1y a2y b2y c2y dyy];
%     
%     N = 100;                % Number of points in x and y
%     Lx = 1e-3;      Ly = Lx;
%     dx = Lx/(N-1);  dy = dx;
%     T_an_1Dprod = zeros(N,N);
% 
%     cd(plot_tool_dir)
%     for i = 1:N
%         x = i*dx;
%         for j = 1:N
%             y = j*dy;
%             T_an_1Dprod(j,i) = Gauss2by1DProduct_discrete(px, py, P0, theta, x, y);
%         end
%     end
%     cd(curr_dir)
%     
%     figure(5)
%     X = 0:dx:Lx;
%     Y = 0:dy:Ly;
% %     contourf(X, Y, T_an_1Dprod, [350 400 500 600 700 1000 1200 1500 1800 1923 3315 5000 8000 10000]);
%     contourf(X, Y, T_an_1Dprod, 14);
%     colorbar
% %     caxis([350, 3500])
% %     caxis manual;
%     title(['t = ', num2str(t-1)]);

    
    %----- 1D 2-term Plotting in 2D -----
    d1d = 0.5*(dxx1d+dyy1d);
    p = [a1x1d b1x1d c1x1d a2x1d b2x1d c2x1d a1y1d b1y1d c1y1d a2y1d b2y1d c2y1d d1d];
    
    % Calculate T(x,y) based on analytical function
    Np = 100;                % Number of points in x and y
    Lx = 1e-3;      Ly = Lx;
    dx = Lx/(Np-1);  dy = dx;
    T_an = zeros(Np,Np);
    X = 0:dx:Lx;
    Y = 0:dy:Ly;
    cd(plot_tool_dir)
    for i = 1:Np
        x = i*dx;
        for j = 1:Np
            y = j*dy;
            T_an(j,i) = Gauss2Eq2DDiscreteFunc(p, P0, theta, x, y);
        end
    end
    cd(curr_dir)
    
    % Plot 2D contour plot of temperature
    figure(3)
    subplot(1,3,2)
    %     colormap(jet);
    contourf(X, Y, T_an, [350 400 500 600 700 1000 1200 1500 1800 1923 3315 5000 8000 10000]);
    colorbar
    caxis([300, 3500])
    caxis manual;
    title(['t = ', num2str(t-1)]);
    xlabel('X'), ylabel('Y')
    title('T_f_i_t_,_1_D')
    
    
    %----- 2D Fit Params -----     
    filename3 = ['Result_T', num2str(t-1), '_Gauss_2DParams.csv'];
    N = csvread(filename3, 1);

    a1x = N(:,1);        b1x = N(:,2);        c1x = N(:,3);
    a2x = N(:,4);        b2x = N(:,5);        c2x = N(:,6);
    
    a1y = N(:,7);        b1y = N(:,8);        c1y = N(:,9);
    a2y = N(:,10);       b2y = N(:,11);       c2y = N(:,12);
    d = N(:,13);
    P0(1) = N(:,14);    P0(2) = N(:,15);    theta = N(:,16);
    
    p = [a1x b1x c1x a2x b2x c2x a1y b1y c1y a2y b2y c2y d];
        cd(plot_tool_dir)
    for i = 1:Np
        x = i*dx;
        for j = 1:Np
            y = j*dy;
            T_an(j,i) = Gauss2Eq2DDiscreteFunc(p, P0, theta, x, y);
        end
    end
    cd(curr_dir)
    
    figure(3)
    subplot(1,3,3)
    %     colormap(jet);
    contourf(X, Y, T_an, [350 400 500 600 700 1000 1200 1500 1800 1923 3315 5000 8000 10000]);
    colorbar
    caxis([300, 3500])
    caxis manual;
    title(['t = ', num2str(t-1)]);
    xlabel('X'), ylabel('Y')
    title('T_f_i_t_,_2_D')
    
    %----- SallyReader- Actual temperature plots -----
    figure(3)
    subplot(1,3,1)
    SallyT=csvread(csvfiles(t).name);
    p=SallyT(:,1:2);
    T2=SallyT(:,4);
    F = TriScatteredInterp(p,T2);
    qz = F(qx,qy);
    %             mesh(qx,qy,qz);
    
    colormap(jet);
    contourf(qx,qy,qz, [350 400 500 600 700 1000 1200 1500 1800 1923 3315 5000 8000 10000]), hold on;
    caxis([300, 3500])
    caxis manual;
    xlabel('X'), ylabel('Y')
    title('T_s_i_m')

    %             plot(SallyT(:,1), SallyT(:,2), '*')
    hold off
    %              caxis([(300)/max(max(T)) (1500)/max(max(T))]*max(max(T)));
    
    colorbar
%     title('Thermal contour')
%     xlabel('Top Surface Domain  in the x direction');
%     ylabel('Top Surface Domain in the y direction')
    %             pause(0.1)
    %     plot(SallyT(:,1), SallyT(:,2), '*'), hold off
    pause
%     hold off
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
xlabel('x'' [m]'), ylabel('T [K]');
hold off

subplot(1,2,2)
plot(M(x_end+1:end,2), M(x_end+1:end,4), 'ro')
axis([-4e-4 4e-4 0 7e4])
xlabel('y'' [m]'), ylabel('T [K]');
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

