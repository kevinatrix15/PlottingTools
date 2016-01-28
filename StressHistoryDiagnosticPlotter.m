%% Read in stress history diagnostic points, and plot data
clear all, clc

data_dir = pwd;
plot_tool_dir = 'C:\Users\Kevin\OneDrive\3DSIM';
cd(data_dir);

filename = ['Result_L1_StressHistoryDiagnostic.csv'];
N = csvread(filename, 1);

xlocs = N(:,1);
ylocs = N(:,2);
zlocs = N(:,3);
times = N(:,4);
sigxx = N(:,5);
sigyy = N(:,6);
sigxy = N(:,7);

num_pts = 110;
num_times = length(times)/num_pts;



sigxx_by_point = zeros(num_times, num_pts);
for t=1:num_times
    sigxx_by_point(t,:) = sigxx((t-1)*num_pts + 1:(t)*num_pts);
end

Lx = max(xlocs)
dx = Lx/num_pts
Ly = max(ylocs)
dy = Ly/num_pts

% XX = 1:num_pts
d_dist = sqrt(dx^2 + dy^2);
XX = 0.5*d_dist:d_dist:sqrt(Lx^2+Ly^2)-0.5*d_dist

figure(1)
plot(XX, sigxx_by_point(146,:), 'ob')
xlabel('Dist. from lower rt. corner [m]', 'FontSize', 14)
ylabel('$$\sigma_{xx}$$ [Pa]','Interpreter','latex', 'FontSize', 14)


figure(2)
plot(times(1:num_pts:end), sigxx_by_point(:,66))
dt = 5;
% step = 567/5

%% Read in stress history diagnostic points, and plot data- both diags
clear all, clc

data_dir = pwd;
plot_tool_dir = 'C:\Users\Kevin\OneDrive\3DSIM';
cd(data_dir);

filename = ['Result_L1_StressHistoryDiagnostic.csv'];
N = csvread(filename, 1);

% 1 suffix is for diagonal from bottom left-top right, 2- bottom right to top left
xlocs1 = N(1:2:end,1);          
xlocs2 = N(2:2:end,1);
ylocs = N(1:2:end,2);
zlocs = N(1:2:end,3);
times = N(1:2:end,4);
sigxx1 = N(1:2:end,5);
sigxx2 = N(2:2:end,5);
sigyy1 = N(1:2:end,6);
sigyy2 = N(2:2:end,6);
sigxy1 = N(1:2:end,7);
sigxy2 = N(2:2:end,7);

num_pts = 110;
num_times = length(times)/num_pts;

XX = 1:num_pts

sigxx1_by_point = zeros(num_times, num_pts);
for t=1:num_times
    sigxx1_by_point(t,:) = sigxx1((t-1)*num_pts + 1:(t)*num_pts);
end

Lx = max(xlocs1)
dx = Lx/num_pts
Ly = max(ylocs)
dy = Ly/num_pts

figure(1)
plot(XX, sigxx1_by_point(146,:))


figure(2)
plot(times(1:num_pts:end), sigxx_by_point(:,66))
dt = 5;
