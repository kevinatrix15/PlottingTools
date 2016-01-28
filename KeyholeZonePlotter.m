%% This script is used to plot the 'mid plane' of a melt pool, showing
%% locations of front and back zones for comparison with matlab code
clear all, clc, close all
tot_steps = 600;

plot_tool_dir = 'C:\Users\Kevin\OneDrive\3DSIM'
curr_dir = pwd; 

for t=534:tot_steps+1
    filename = ['Result_L1_O16_T', num2str(t-1), '_Keyhole.csv'];
   
    cd(plot_tool_dir)
    [x_liq, z_liq, x_vap, z_vap, x_lqcor, z_lqcor, x_vapcor, z_vapcor,...
        x_bkvel, z_bkvel, u_bkvel, x_frvel, z_frvel, u_frvel] = readKeyholeFiles(filename, curr_dir);
%     [x_vap, z_vap] = readKeyholeLocations(filename2, curr_dir);
%     [x_fr, z_fr] = readKeyholeLocations(filename3, curr_dir);
%     [x_bk, z_bk] = readKeyholeLocations(filename4, curr_dir);
% %     [x_frZ, z_frZ] = readKeyholeLocations(filename5, curr_dir);
% %     [x_bkZ, z_bkZ] = readKeyholeLocations(filename6, curr_dir);
%     [x_lqcor, z_lqcor] = readKeyholeLocations(liqCornerFile, curr_dir);
%     [x_vapcor, z_vapcor] = readKeyholeLocations(vapCornerFile, curr_dir);
%     [x_frvel, z_frvel, u_frvel] = readKeyholeVels(velfile2, curr_dir);
%     [x_bkvel, z_bkvel, u_bkvel] = readKeyholeVels(velfile1, curr_dir);
%     cd(curr_dir)
    
    x_lqcor(end+1) = x_lqcor(1);
    z_lqcor(end+1) = z_lqcor(1);
    x_vapcor(end+1) = x_vapcor(1);
    z_vapcor(end+1) = z_vapcor(1);
    
    figure(1)
    plot(x_liq, z_liq, '-bo', x_vap, z_vap, '-ro', 'LineWidth', 3)
    hold on
    plot(x_lqcor, z_lqcor, '--cx', x_vapcor, z_vapcor, '--cx', 'LineWidth', 2, 'MarkerSize', 18), hold on
    plot(x_vapcor, z_vapcor, '--yx', 'LineWidth', 2, 'MarkerSize', 18), hold on   
%     plot(x_frZ, z_frZ, 'gx', x_bkZ, z_bkZ, 'mx', 'LineWidth', 2)
    quiver(x_bkvel, z_bkvel, u_bkvel, zeros(length(u_bkvel),1))
    quiver(x_frvel, z_frvel, u_frvel, zeros(length(u_frvel),1))
    hold off
    axis([1e-4 8e-4 4.8e-3 5.04e-3])
    xlabel('X [m]'), ylabel('Z [m]')
    title(['t_s_t_e_p = ',num2str(t)])
    pause
end

%%
clear all, clc, close all
tot_steps = 900;

plot_tool_dir = 'C:\Users\Kevin\OneDrive\3DSIM'
curr_dir = pwd; 

for t=1:tot_steps+1
    filename1 = ['Result_L1_T', num2str(t-1), '_LiqBoundary.csv'];
    filename2 = ['Result_L1_T', num2str(t-1), '_VapBoundary.csv'];
    filename3 = ['Result_L1_T', num2str(t-1), '_FrontHull.csv'];
    filename4 = ['Result_L1_T', num2str(t-1), '_BackHull.csv'];
%     filename5 = ['Result_L1_T', num2str(t-1), '_FrontZonePts.csv'];
%     filename6 = ['Result_L1_T', num2str(t-1), '_BackZonePts.csv'];
    velfile1 = ['Result_L1_T', num2str(t-1), '_BackVels.csv'];
    velfile2 = ['Result_L1_T', num2str(t-1), '_FrontVels.csv'];
    vapCornerFile = ['Result_L1_T', num2str(t-1), '_VapCorners.csv'];
    liqCornerFile = ['Result_L1_T', num2str(t-1), '_LiqCorners.csv'];
    
    cd(plot_tool_dir)
    [x_liq, z_liq] = readKeyholeLocations(filename1, curr_dir);
    [x_vap, z_vap] = readKeyholeLocations(filename2, curr_dir);
    [x_fr, z_fr] = readKeyholeLocations(filename3, curr_dir);
    [x_bk, z_bk] = readKeyholeLocations(filename4, curr_dir);
%     [x_frZ, z_frZ] = readKeyholeLocations(filename5, curr_dir);
%     [x_bkZ, z_bkZ] = readKeyholeLocations(filename6, curr_dir);
    [x_lqcor, z_lqcor] = readKeyholeLocations(liqCornerFile, curr_dir);
    [x_vapcor, z_vapcor] = readKeyholeLocations(vapCornerFile, curr_dir);
    [x_frvel, z_frvel, u_frvel] = readKeyholeVels(velfile2, curr_dir);
    [x_bkvel, z_bkvel, u_bkvel] = readKeyholeVels(velfile1, curr_dir);
    cd(curr_dir)
    
    x_lqcor(end+1) = x_lqcor(1);
    z_lqcor(end+1) = z_lqcor(1);
    x_vapcor(end+1) = x_vapcor(1);
    z_vapcor(end+1) = z_vapcor(1);
    
    figure(1)
    plot(x_liq, z_liq, '-bo', x_vap, z_vap, '-ro', 'LineWidth', 3)
    hold on
    plot(x_lqcor, z_lqcor, '--yx', x_vapcor, z_vapcor, '--cx', 'LineWidth', 2, 'MarkerSize', 18), hold on
    plot(x_vapcor, z_vapcor, '--cx', 'LineWidth', 2, 'MarkerSize', 18), hold on   
    plot(x_fr, z_fr, '-gx', x_bk, z_bk, '-mx', 'LineWidth', 2)
%     plot(x_frZ, z_frZ, 'gx', x_bkZ, z_bkZ, 'mx', 'LineWidth', 2)
    quiver(x_bkvel, z_bkvel, u_bkvel, zeros(length(u_bkvel),1))
    quiver(x_frvel, z_frvel, u_frvel, zeros(length(u_frvel),1))
    hold off
    axis([1e-4 8e-4 4.8e-3 5.04e-3])
    xlabel('X [m]'), ylabel('Z [m]')
    title(['t_s_t_e_p = ',num2str(t)])
    pause
end
