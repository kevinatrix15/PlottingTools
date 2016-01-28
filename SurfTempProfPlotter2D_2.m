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

fitting_method = 'SkewNorm';    % 'Gauss2', 'Gauss1'

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
%     axis([-6e-4 6e-4 0 5000])
%     axis([-4e-4 4e-4 0 4000])
    xlabel('x'' [m]'), ylabel('T [K]');
    title(['t = ', num2str(t-1)])
    
    
    yy = M(x_end+1:end,2), yprof = M(x_end+1:end,4);
    subplot(1,2,2)
    plot(yy, yprof, 'ro'), hold on
%     axis([-6e-4 6e-4 0 5000])
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
%     axis([-6e-4 6e-4 0 5000])
    %     pause(0.5)
    %     hold off
    %     figure(2)
    %     plot(M(x_end+1:end,2), M(x_end+1:end,4), '-x')
    
    if (strcmp(fitting_method, 'Gauss2'))
        filename2 = ['Result_T', num2str(t-1), '_Gauss2_1DParams.csv'];
        N = csvread(filename2, 1);
        a1x1d = N(:,1);        b1x1d = N(:,2);        c1x1d = N(:,3);
        a2x1d = N(:,4);        b2x1d = N(:,5);        c2x1d = N(:,6);
        dxx1d = N(:,7);
        
        a1y1d = N(:,8);        b1y1d = N(:,9);        c1y1d = N(:,10);
        a2y1d = N(:,11);       b2y1d = N(:,12);       c2y1d = N(:,13);
        dyy1d = N(:,14);
        
        %     P0(1) = N(:,15);    P0(2) = N(:,16);    theta = N(:,17);
        
        figure(2)
        subplot(1,2,1)
        Gauss2Eq1DPlotter(a1x1d, b1x1d, c1x1d, a2x1d, b2x1d, c2x1d, dxx1d, -6e-4, 6e-4,'--k');
        
        subplot(1,2,2)
        Gauss2Eq1DPlotter(a1y1d, b1y1d, c1y1d, a2y1d, b2y1d, c2y1d, dyy1d, -6e-4, 6e-4, '--k');
        %     axis([-3e-4 3e-4 0 5000])
        %pause(0.2)
        hold off
        
    elseif (strcmp(fitting_method, 'SkewNorm'))
        filename2 = ['Result_T', num2str(t-1), '_SkewNorm2_1DParams.csv'];
        N = csvread(filename2, 1);
        ax = N(:,1);    bx = N(:,2);    cx = N(:,3);    ddx = N(:,4);    alphx = N(:,5);
        ay = N(:,6);    by = N(:,7);    cy = N(:,8);    ddy = N(:,9);    alphy = N(:,10);
        P0(1) = N(:,11);    P0(2) = N(:,12);    theta = N(:,13);

        %     P0(1) = N(:,15);    P0(2) = N(:,16);    theta = N(:,17);
        
        figure(2)
        subplot(1,2,1)
        SkewNorm1DPlotter(ax, bx, cx, ddx, alphx, -6e-4, 6e-4,'--k');
        hold off
        
        subplot(1,2,2)
        SkewNorm1DPlotter(ay, by, cy, ddy, alphy, -6e-4, 6e-4,'--k');
        hold off
        pause
        
        continue
    end
%     axis([-4e-4 2e-4 0 5000])
    legend('T_s_i_m', 'T_f_i_t')
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

    

    %----- 1D coefficient Plotting in 2D -----
    if (strcmp(fitting_method, 'Gauss2'))
        filename3 = ['Result_T', num2str(t-1), '_Gauss_2DParams.csv'];
        N = csvread(filename3, 1);
        P0(1) = N(:,14);    P0(2) = N(:,15);    theta = N(:,16);
        
        d1d = 0.5*(dxx1d+dyy1d);
        p = [a1x1d/2 b1x1d c1x1d a2x1d/2 b2x1d c2x1d a1y1d/2 b1y1d c1y1d a2y1d/2 b2y1d c2y1d d1d];
        
    elseif (strcmp(fitting_method, 'SkewNorm'))
        dmean = 0.5*(ddx + ddy);
        p = [ax/2 bx cx alphx ay/2 by cy alphy dmean];
    end
    
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
            if (strcmp(fitting_method, 'Gauss2'))
                T_an(j,i) = Gauss2Eq2DDiscreteFunc(p, P0, theta, x, y);
            elseif (strcmp(fitting_method, 'SkewNorm'))
                T_an(j,i) = Gauss2Eq2DDiscreteFunc(p, P0, theta, x, y);
            end
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
    
    
%     figure(11)
%     plot(b1x1d,b1y1d,'*r',b1x1d,b2y1d,'ob',b2x1d,b1y1d,'kd', b2x1d,b2y1d, 'sm')
%     legend('b1x,b1y', 'b1x,b2y', 'b2x,b1y', 'b2x,b2y')
%     hold on
%     
%     ellipse(b1x1d, b1y1d, c1x1d, c1y1d, 'r')
%     
%     ellipse(b1x1d, b2y1d, c1x1d, c2y1d, 'b')
%     
%     ellipse(b2x1d, b1y1d, c2x1d, c1y1d, 'k')
%     
%     ellipse(b2x1d, b2y1d, c2x1d, c2y1d, 'm')
%         
%     axis([-1e-4 1e-4 -1e-4 1e-4])
%     asum = a1x1d+a2x1d+a1y1d+a2y1d;
%     text(5e-5, 5e-5, ['a1x = ', num2str(a1x1d)]);
%     text(5e-5, 4e-5, ['a2y = ', num2str(a2y1d)]);
%     text(5e-5, 3e-5, ['a1y = ', num2str(a1y1d)]);
%     text(5e-5, 2e-5, ['a2x = ', num2str(a2x1d)]);
%     text(5e-5, 1e-5, ['Sum(a) = ', num2str(asum)]);
% 
%     hold off

    %----- Coefficient-corrected 2D function (using 1D params) -----%
    p = [a1x1d b1x1d c1x1d a2x1d b2x1d c2x1d a1y1d b1y1d c1y1d a2y1d b2y1d c2y1d d1d];
    pCC = fittingCoeffSwap(p);      % Returns corrected coefficient values for 2D plot

    T_an_CC = zeros(Np,Np);
    X = 0:dx:Lx;
    Y = 0:dy:Ly;
    cd(plot_tool_dir)
    for i = 1:Np
        x = i*dx;
        for j = 1:Np
            y = j*dy;
            T_an_CC(j,i) = genCoeffGauss2Eq2D(pCC, P0, theta, x, y);
        end
    end
    cd(curr_dir)
    
    % Plot 2D contour plot of temperature
    figure(3)
    subplot(1,3,3)
    %     colormap(jet);
    contourf(X, Y, T_an_CC, [350 400 500 600 700 1000 1200 1500 1800 1923 3315 5000 8000 10000]);
    colorbar
    caxis([300, 3500])
    caxis manual;
    title(['t = ', num2str(t-1)]);
    xlabel('X'), ylabel('Y')
    title('T_f_i_t_,_C_C')
    



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
    
    
%     %----- Matlab Levenberg-Marquardt fit plot -----
%     filename4 = ['Result_T', num2str(t-1), '_InterpolatedTemps2D.csv'];
%     MM = csvread(filename4, 1);
%     
%     xyLocs(:,1) = MM(:,1); xyLocs(:,2) = MM(:,2); Tvals = MM(:,3);
%     
%     plot_tool_dir = 'C:\Users\Kevin\OneDrive\3DSIM'
%     cd(plot_tool_dir)
%     
% %     pLM_0 = [a1x1d b1x1d c1x1d a2x1d b2x1d c2x1d a1y1d b1y1d c1y1d a2y1d b2y1d c2y1d d1d];
%     pLM_0 = [max(xprof) 1e-6 1e-6 max(xprof) 1e-6 1e-6 max(xprof) 1e-6 1e-6 max(xprof) 1e-6 1e-6 353];
% %     x0 = [max(xprof); 0; 3e-6; 3733.32; -0.00020238; 1e-6; min(xprof)];
%     options = optimset('Algorithm','levenberg-marquardt');
%     [pLM, resnorm, residual, exitflag, output] = lsqcurvefit(@Gauss2Eq2D_funFit, pLM_0, xyLocs, Tvals, [], [], options);
%     
%     % Plot Matlab Fit
%     T_an_LM = zeros(Np,Np);
%     cd(plot_tool_dir)
%     for i = 1:Np
%         x = i*dx;
%         for j = 1:Np
%             y = j*dy;
%             T_an_LM(j,i) = Gauss2Eq2DDiscreteFunc(pLM, P0, theta, x, y);
%         end
%     end
%     cd(curr_dir)
%     
%     figure(3)
%     subplot(1,3,2)
%     %     colormap(jet);
%     contourf(X, Y, T_an_LM, [350 400 500 600 700 1000 1200 1500 1800 1923 3315 5000 8000 10000]);
%     colorbar
%     caxis([300, 3500])
%     caxis manual;
%     title(['t = ', num2str(t-1)]);
%     xlabel('X'), ylabel('Y')
%     title('T_f_i_t_,_L_M')
    
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
%     pause(0.1)
    
    figure(6)
    plot(t,a1x,'*r',t,a2x,'ob',t,a1y,'^k',t,a2y,'sm')
    legend('a1x','a2x','a1y','a2y')
    hold on
    
    figure(7)
    subplot(1,2,1)
    plot(t,a1x,'*r',t,a2x,'ob',t,a1x+a2x,'kd')
    legend('a1x','a2x','a1x+a2x')
    hold on
     
    subplot(1,2,2)
    plot(t,a1y,'*r',t,a2y,'ob',t,a1y+a2y,'kd')
    legend('a1y','a2y','a1y+a2y')
    hold on
    
    figure(8)
    plot(b1x,b1y,'*r',b1x,b2y,'ob',b2x,b1y,'kd', b2x,b2y, 'sm')
    legend('b1x,b1y', 'b1x,b2y', 'b2x,b1y', 'b2x,b2y')
    hold on

    figure(9)
    xSq = [b1x, b2x, b2x, b1x, b1x];
    ySq = [b1y, b1y, b2y, b2y, b1y];
    plot(xSq, ySq, '-^k')
    legend('b1x,b1y', 'b1x,b2y', 'b2x,b1y', 'b2x,b2y')
    hold on
%     hold off
    figure(10)
    subplot(1,2,1)
    plot(t,b1x,'*r',t,b2x,'ob')
    legend('b1x','b2x')
    hold on
     
    subplot(1,2,2)
    plot(t,b1y,'*r',t,b2y,'ob')
    legend('b1y','b2y')
    hold on
    
    
%     figure(11)
%     plot(b1x,b1y,'*r',b1x,b2y,'ob',b2x,b1y,'kd', b2x,b2y, 'sm')
%     legend('b1x,b1y', 'b1x,b2y', 'b2x,b1y', 'b2x,b2y')
%     hold on
%     
%     ellipse(b1x, b1y, c1x, c1y, 'r')
%     
%     ellipse(b1x, b2y, c1x, c2y, 'b')
%     
%     ellipse(b2x, b1y, c2x, c1y, 'k')
%     
%     ellipse(b2x, b2y, c2x, c2y, 'm')
%         
%     axis([-1e-4 1e-4 -1e-4 1e-4])
%     asum = a1x+a2x+a1y+a2y;
%     text(5e-5, 5e-5, ['a1x = ', num2str(a1x)]);
%     text(5e-5, 4e-5, ['a2y = ', num2str(a2y)]);
%     text(5e-5, 3e-5, ['a1y = ', num2str(a1y)]);
%     text(5e-5, 2e-5, ['a2x = ', num2str(a2x)]);
%     text(5e-5, 1e-5, ['Sum(a) = ', num2str(asum)]);
% 
%     hold off
    
    figure(12)
    subplot(1,2,1)
    plot(t,c1x1d,'*r',t,c2x1d,'ob')
    legend('c1x','c2x')
    hold on
     
    subplot(1,2,2)
    plot(t,c1y1d,'*r',t,c2y1d,'ob')
    legend('c1y','c2y')
    hold on
    
    
    figure(13)
    subplot(1,2,1)
    plot(a1x1d, abs(b2x1d-b1x1d)/(c1x1d+c2x1d), '*r')
    ylabel('|b2x - b1x|/(c1x+c2x)'), xlabel('a1x')
    hold on
    
    subplot(1,2,2)
    plot(a1y1d, abs(b2y1d-b1y1d)/(c1y1d+c2y1d), '*r')
    ylabel('|b2y - b1y|/(c1y+c2y)'), xlabel('a1y')
    hold on
    pause
end


%% This cell reads in 2D interpolated (simulation) temperature values at
% top surface and performs nonlinear LSQ fit to data

clear all, clc, close all
tot_steps = 1800;
curr_dir = pwd;

for t=1:1:tot_steps+1
    
    %----- 1D 2 term Ploting -----
    filename4 = ['Result_T', num2str(t-1), '_InterpolatedTemps2D.csv'];
    M = csvread(filename4, 1);
    
    xLocs = M(:,1); yLocs = M(:,2); Tvals = M(:,3);
    
    plot_tool_dir = 'C:\Users\Kevin\OneDrive\3DSIM'
    cd(plot_tool_dir)
    pLM = [a1x1d b1x1d c1x1d a2x1d b2x1d c2x1d a1y1d b1y1d c1y1d a2y1d b2y1d c2y1d d1d];
%     x0 = [max(xprof); 0; 3e-6; 3733.32; -0.00020238; 1e-6; min(xprof)];
    options = optimset('Algorithm','levenberg-marquardt');
    [X, resnorm, residual, exitflag, output] = lsqcurvefit(@Gauss2Eq2D_funFit, pLM, xLocs, yLocs, Tvals, [], [], options);
    figure(20)
    subplot(1,2,1)
%     Gauss2Eq1DPlotter(X(1), X(2), X(3), X(4), X(5), X(6), X(7), -4e-4, 4e-4, '-g')
    cd(curr_dir)
    

    xx = M(1:x_end,1), xprof = M(1:x_end,4);
    figure(2)
    subplot(1,2,1)
    plot(xx, xprof, 'bo'), hold on
%     axis([-6e-4 6e-4 0 5000])
%     axis([-4e-4 4e-4 0 4000])
    xlabel('x'' [m]'), ylabel('T [K]');
    title(['t = ', num2str(t-1)])
    
    
    yy = M(x_end+1:end,2), yprof = M(x_end+1:end,4);
    subplot(1,2,2)
    plot(yy, yprof, 'ro'), hold on
%     axis([-6e-4 6e-4 0 5000])
    %axis([-4e-4 4e-4 0 4000])
    xlabel('y'' [m]'), ylabel('T [K]');
end

    
    
