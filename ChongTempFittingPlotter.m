%% This script parses and plots surface temperatre profiles along scan-aligned axis directions (x', y')

% Chong- set directory to the simulation Output directory

clear all, clc, close all
tot_steps = 1800;
curr_dir = pwd;

%----- The below is for SallyReader -----
Nodes=26;
sd=8;
SallyT=csvread('Result_L1_O1_T0_Temperature.csv');
bedwidth_x=max(SallyT(:,1));
bedwidth_y=max(SallyT(:,2));
[qx,qy] = meshgrid([0:bedwidth_x/(sd*(Nodes-1)):bedwidth_x],[0:bedwidth_y/(sd*(Nodes-1)):bedwidth_y]);
csvfiles = dir('*Temperature.csv'); 
[~,index] = sortrows([csvfiles.datenum].'); csvfiles = csvfiles(index);
numfiles = length(csvfiles);

for t=1:10:tot_steps+1
    
    %----- Interpolated Temperature along x' and y' -----
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
    xlabel('x'' [m]'), ylabel('T [C]');
    title(['t = ', num2str(t-1)])
    
    
    yy = M(x_end+1:end,2), yprof = M(x_end+1:end,4);
    subplot(1,2,2)
    plot(yy, yprof, 'ro'), hold on
    axis([-6e-4 6e-4 0 5000])
    %axis([-4e-4 4e-4 0 4000])
    xlabel('y'' [m]'), ylabel('T [C]');
    
    
    %----- 1D fit for x' and y' -----
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
    Gauss2Eq1DPlotter(a1x, b1x, c1x, a2x, b2x, c2x, dx, -6e-4, 6e-4,'--k');
    axis([-6e-4 6e-4 0 5000])
    hold off
    
    subplot(1,2,2)
    Gauss2Eq1DPlotter(a1y, b1y, c1y, a2y, b2y, c2y, dy, -6e-4, 6e-4, '--k');
    axis([-6e-4 6e-4 0 5000])
    %pause(0.2)
    hold off
    
    
    %----- Chong- TO DO: sort peaks according to method we discussed, and
    % use/modify Gauss2Eq2DDiscreteFunc.m to create 2D fit -----
    
    %                       Do stuff here...
    
    %----------------------------------------------------------------------
    
    
    %----- SallyReader- Actual 2D temperature plots -----
    figure(3)
%     subplot(1,2,2)
    SallyT=csvread(csvfiles(t).name);
    p=SallyT(:,1:2);
    T2=SallyT(:,4);
    F = TriScatteredInterp(p,T2);
    qz = F(qx,qy);
    %             mesh(qx,qy,qz);
    
    colormap(jet);
    contourf(qx,qy,qz, [350 400 500 600 700 1000 1200 1500 1800 1923 3315 5000 8000 10000]), hold on;
    caxis([350, 3500])
    caxis manual;
    %             plot(SallyT(:,1), SallyT(:,2), '*')
    hold off
    %              caxis([(300)/max(max(T)) (1500)/max(max(T))]*max(max(T)));
    
    colorbar
    title('Thermal contour')
    xlabel('Top Surface Domain  in the x direction');
    ylabel('Top Surface Domain in the y direction')
    %             pause(0.1)
    %     plot(SallyT(:,1), SallyT(:,2), '*'), hold off
    pause(0.25)
%     hold off
end