%% Plot cumulative displacements due to thermal expansion / contraction
clear all, clc

data_dir = pwd;
plot_tool_dir = 'C:\Users\Kevin\OneDrive\3DSIM';
cd(data_dir);

num_layers = 1;             % Number of layers considered by diagnostic
num_times = 30000;

for l = 1:num_layers
    for tt = 1:1:3600
    filename = ['Result_L', num2str(l), '_T', num2str(tt), '_ShrinkageDisplacements.csv'];

    N = csvread(filename, 1);

    xlocs = N(:,1);
    ylocs = N(:,2);
    xdisp = N(:,3);
    ydisp = N(:,4);
    
    num_pts = length(xlocs);
    
    nx = find(ylocs(2:end) ~= ylocs(1), 1, 'first');
    ny = num_pts / nx;

    dx = xlocs(2) - xlocs(1);
    dy = ylocs(nx+1) - ylocs(1);

    Lx = xlocs(end);
    Ly = ylocs(end);

    X = 0:dx:Lx;
    Y = 0:dy:Ly;
    [YY,XX] = meshgrid(Y,X);

    U = zeros(nx,ny);
    V = zeros(nx,ny);
    for i= 1:nx
        for j=1:ny
              U(i,j) = xdisp(i + (j-1)*nx);
              V(i,j) = ydisp(i + (j-1)*nx);
        end
    end

    figure(1)
    % contourf(XX, YY, T_norm(:,:,2)), shading flat, colorbar,colormap('cool'), hold on
    % caxis([0 0.12]);
    quiver( XX(1:1:end,1:1:end), YY(1:1:end,1:1:end), U(1:1:end,1:1:end,1), ...
        V(1:1:end,1:1:end,1), 3, 'k' )
    xlabel('x/H'), ylabel('y/H')
    title(['time =', num2str(tt)])
    
%     pause(0.1)
%     pause
    end
end

