%% Phase Change Diagnostic plotter
clear all, close all, clc

data_dir = pwd;
plot_tool_dir = 'C:\Users\Kevin\OneDrive\3DSIM';
cd(data_dir);

num_layers = 3;             % Number of layers considered by diagnostic

cooling_time = 0.01;        % Cooling time between layers (s)
T_liq = 1690.15;            % Liquidus temperature (K)
T_solidus = 1567.15;        % Solidus temperature (K)
T_ssTrans = 1301.15;        % Solid state transition temperature (K)

for l = 1:num_layers
%     filename = ['Result_L', num2str(l), '_PointTemp_dTdt.csv'];
    filename = ['GE_SolverResult_335_1800.csv_L', num2str(l), '_PointTemp_dTdt.csv'];

    
    % filename = ['Result_L1', '_PointTemp_dTdt.csv'];
    % filename = ['Result_L1', '_StressHistoryDiagnostic.csv'];
        
    N = csvread(filename, 1);
    
    if (l==1)
        xlocs = N(:,1);
        ylocs = N(:,2);
        zlocs = N(:,3);
        times = N(:,4);
        T = N(:,5);
        dTdt = N(:,6);
        
        num_pts = find(times(2:end) ~= times(1), 1, 'first');
        
        nx = find(ylocs(2:end) ~= ylocs(1), 1, 'first');

    else
        times = [times; N(:,4)+times(end)+cooling_time];
        T = [T; N(:,5)];
        dTdt = [dTdt; N(:,6)];
    end
    
    num_times = length(times)/num_pts;
    T_by_point = zeros(num_times, num_pts);
    dTdt_by_point = zeros(num_times, num_pts);
    
    % store time at end of layer
    t_end_layer(l) = times(end);
end

for i=1:num_pts
    for j = 1:num_times
        row = i + (j-1)*num_pts;
        T_by_point(j, i) = T(row);
        dTdt_by_point(j, i) = dTdt(row);
    end
end

%----- Plots -----
% Define Points of interest
is = [35 75 75 35 55];
js = [35 35 75 75 55];
for i = 1:length(is)
    idx(i) = is(i) + (js(i)-1)*nx;
end

% figure(1)
% plot(times(1:num_pts:end), T_by_point, '-o'), hold on
% % plot(times(1:num_pts:end), dTdt_by_point/1e4, '--')
% legend('point 1', 'point 2', 'point 3', 'point 4')
% 
% figure(2)
% plot(times(1:num_pts:end), dTdt_by_point)
% legend('point 1', 'point 2', 'point 3', 'point 4')

figure(3)
lstyles = {'-r', '-g', '-b', '-m', '-k'};
for i=1:length(idx)
    plot(times(1:num_pts:end), T_by_point(:,idx(i)), lstyles{i}, 'LineWidth', 2), hold on
    legend('point 1', 'point 2', 'point 3', 'point 4', 'point 5')
end
plot([0 times(end)], [T_ssTrans T_ssTrans], [0 times(end)], [T_solidus T_solidus], [0 times(end)], [T_liq T_liq]) 

% Vertical lines marking cooling times
for l=1:num_layers-1
    plot([t_end_layer(l) t_end_layer(l)], [0 max(T)], '--k',...
        [t_end_layer(l)+cooling_time t_end_layer(l)+cooling_time],...
        [0 max(T)], '--k') 
end



figure(4)
plot(times(1:num_pts:end), dTdt_by_point(:,5975), '-ro'), hold on
plot(times(1:num_pts:end), dTdt_by_point(:,5990), '-go')
plot(times(1:num_pts:end), dTdt_by_point(:,6000), '-bo')
plot(times(1:num_pts:end), dTdt_by_point(:,6015), '-mo')
% plot(times(1:num_pts:end), dTdt_by_point(:,862), 'ko')

%% Read in phase map
close all 

data_dir = pwd;
plot_tool_dir = 'C:\Users\Kevin\OneDrive\3DSIM';
cd(data_dir);

% filename = ['Result_L3', '_PhaseMap.csv'];
filename = ['GE_SolverResult_335_1800.csv_L1_PhaseMap.csv'];


fid = fopen(filename, 'rt')

header = fgetl(fid);

% M = textscan(fid, '%f %f %f %f %d %f %d %f %d %f %d %f %d', 'delimiter', ',');
M = csvread(filename, 1);

X = M(:,1);
Y = M(:,2);
Z = M(:,3);
% num_pts = find(times(2:end) ~= times(1), 1, 'first');
nx = find(Y(2:end) ~= Y(1), 1, 'first');
ny = length(Y)/nx;
% nz = find(Z(2:end) ~= Z(1), 1, 'first')+1
nz = 1;

dx = X(2) - X(1);
dy = Y(nx+1) - Y(1);

XX = X(1:nx);
YY = Y(1:nx:end);

dt = 25e-6;                  % TODO: Define generally...
num_time_steps = 500;       % TODO: Get from c++...

p = struct;

%----- Assign p struct data -----
for j=1:ny
    for i=1:nx
        idx = i + (j-1)*nx;
        p(idx).X = X(idx);
        p(idx).Y = Y(idx);
        p(idx).times = M(idx,4:2:13);
        p(idx).phases = M(idx,5:2:13);
%         for k=1:(size(M,2)-3)
    end
end

%----- Plot phase map for each time step -----
X_vert = zeros(length(X),4);
Y_vert = zeros(length(Y),4);
Z_vert = zeros(length(Z),4);
for i = 1:length(X)
    % Faces on x-y plane
    %     if norm.x(i)==0 && norm.y(i)==0 && (norm.z(i)==1 || norm.z(i)==-1)
    X_vert(i,1) = X(i);
    X_vert(i,2) = X(i)+dx;
    X_vert(i,3) = X(i)+dx;
    X_vert(i,4) = X(i);
    
    Y_vert(i,1) = Y(i);
    Y_vert(i,2) = Y(i);
    Y_vert(i,3) = Y(i)+dy;
    Y_vert(i,4) = Y(i)+dy;
    
    Z_vert(i,:) = Z(i);
end
        
% clr_map = zeros(1,length(X));
% clr_map = zeros(length(X), 3);

%----- Define 'Probe' Points -----
is = [35 75 75 35 55];
js = [35 35 75 75 55];
for i = 1:length(is)
    idx(i) = is(i) + (js(i)-1)*nx;
end

current_phase_map = zeros(ny, nx);
for tidx=1:num_time_steps
    time = (tidx-1)*dt;
    for j=1:ny
        for i=1:nx
            idx = i + (j-1)*nx;
            
            ii = 2;
        while (ii <= size(p(idx).phases, 2) && p(idx).times(ii) <= time ...
            && p(idx).times(ii) > p(idx).times(ii-1))
                current_phase_map(j,i) = p(idx).phases(ii);
%             while (ii <= size(p(idx).phases, 2) && p(idx).times(ii) <= time)
%                 current_phase_map(j,i) = p(idx).phases(ii);
                
                % Set color based on material type
%                 if (p(idx).phases(ii) == 0)     % Powder
%                     clr_map(ii, :) = [0 0 1];
%                 elseif (p(idx).phases(ii) == 2)     % Liquid
%                     clr_map(ii, :) = [1 0 0];
%                 elseif (p(idx).phases(ii) == 1)     % Solid
%                     clr_map(ii, :) = [0 1 0];
% %                 else
% %                     clr_map(ii) = 'w';
%                 end


%                 if (p(idx).phases(ii) == 0)     % Powder
%                     clr_map(ii) = 'b';
%                 elseif (p(idx).phases(ii) == 2)     % Liquid
%                     clr_map(ii) = 'r';
%                 elseif (p(idx).phases(ii) == 1)     % Solid
%                     clr_map(ii) = 'g';
%                 else
%                     clr_map(ii) = 'w';
%                 end
                
                ii = ii+1;
            end
            
        end
    end
%     figure(5)
%     h = fill3(X_vert, Y_vert, Z_vert, clr_map);
%     set(h, 'EdgeAlpha', 0);
    
    figure(6)
    pcolor(current_phase_map);
    colorbar
%     contourf(XX, YY, current_phase_map)
        caxis([0 2])
    caxis manual;
    
    hold on;
    for i=1:length(is)
        plot(is(i), js(i), 'm*', 'MarkerSize', 14, 'LineWidth', 2)
    end
    hold off
    
    title(['Time = ', num2str(time), ' s'])
%     colorbar
    pause(0.1)
end
        
        