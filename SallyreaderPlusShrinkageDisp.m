clear all, clc


%----- Shrinkage Parameters -----
data_dir = pwd;
plot_tool_dir = 'C:\Users\Kevin\OneDrive\3DSIM';
cd(data_dir);

plotShrinkage = true;
% shrinkageFileStem = 'shrinkagefile';
shrinkageFileStem = 'ShrinkageDisplacements';


plotTemperature = false;
temperatureFileStem = 'thermal';
% Tfield_dir = '../Output_102715_Animation';
Tfield_dir = pwd;

num_layers = 2;             % Number of layers considered by diagnostic
num_times = 30000;


N=1000;
windowing=1;
windowing1=20;
Nodes=21;
sd=4;
% mpfile=csvread('Result_L1_Meltpool.csv',1,0);
% max_x=max(mpfile(:,4));
% min_x=min(mpfile(:,4));
% max_y=max(mpfile(:,5));
% min_y=min(mpfile(:,5));
max_x=2.5e-3;
min_x=0;
max_y=2.5e-3;
min_y=0;

bedwidth_x=1.2*(max_x-min_x);
bedwidth_y=1.2*(max_y-min_y);
bedmin_x=min_x;
bedmax_x=max_x+0.2*bedwidth_x;
bedmin_y=min_y;
bedmax_y=max_y+0.2*bedwidth_y;

if (plotTemperature == true)
    cd(Tfield_dir)
    SallyT=csvread([temperatureFileStem, '_L1_O17_T100_Temperature.csv']);
    
    % set(gca, 'nextplot', 'replacechildren');
    [qx,qy] = meshgrid([bedmin_x:bedwidth_x/(sd*(Nodes-1)):bedmax_x],[bedmin_y:bedwidth_y/(sd*(Nodes-1)):bedmax_y]);
    csvfiles = dir('*Temperature.csv');
%     [~,index] = sortrows([csvfiles.datenum].'); csvfiles = csvfiles(index);
    numfiles = length(csvfiles);
    cd(data_dir)
end

saveFigure = false;

% U = zeros(90,90);
% V = zeros(90,90);

for t = 1:1:100
    clf
    ppp = figure(10)
    
    if (plotTemperature == true)
        cd(Tfield_dir)
        SallyT=csvread(csvfiles(t).name);
        cd(data_dir)
        p=SallyT(:,1:2);
        %     plot(p(:,1),p(:,2),'r*')
        T2=SallyT(:,4);
        F = TriScatteredInterp(p,T2);
        qz = F(qx,qy);
        %             mesh(qx,qy,qz);
        
        hold on
        contour(qx,qy,qz, [350 400 500 600 700 1000 1200 1500 1800 1923 3315]);
        caxis([350, 3315])
        caxis manual;
        %              caxis([(300)/max(max(T)) (1500)/max(max(T))]*max(max(T)));
        
        colorbar
        pause
    end
    
    %----- Plot Shrinkage -----
    if plotShrinkage == true
        filename = [shrinkageFileStem, '_L', num2str(t+1000), '.csv'];
        %     filename = ['Result_L', num2str(1), '_T', num2str(t), '_ShrinkageDisplacements.csv'];
        
        cd(plot_tool_dir)
        [U, V, Eps1, X, Y] = readShrinkage(pwd, data_dir, filename);
        cd(data_dir)
        
        %     filename = ['Result_L', num2str(1+1000), '_ShrinkageDisplacements.csv'];
        %     filename = ['Result_L', num2str(1), '_T', num2str(t), '_ShrinkageDisplacements.csv'];
        %     N = csvread(filename, 1);
        %
        %     xlocs = N(:,1);
        %     ylocs = N(:,2);
        %     xdisp = N(:,3);
        %     ydisp = N(:,4);
        %
        %     num_pts = length(xlocs);
        %
        %     nx = find(ylocs(2:end) ~= ylocs(1), 1, 'first');
        %     ny = num_pts / nx;
        %
        %     dx = xlocs(2) - xlocs(1);
        %     dy = ylocs(nx+1) - ylocs(1);
        %
        %     Lx = xlocs(end);
        %     Ly = ylocs(end);
        %
        %     X = 0:dx:Lx;
        %     Y = 0:dy:Ly;
        [YY,XX] = meshgrid(Y,X);
        
        %     U = zeros(nx,ny);
        %     V = zeros(nx,ny);
        %     Uprev = U;
        %     Vprev = V;
        %     for i= 1:nx
        %         for j=1:ny
        %             U(i,j) = xdisp(i + (j-1)*nx);
        %             V(i,j) = ydisp(i + (j-1)*nx);
        %         end
        %     end
        %
        %     figure(1)
        % contourf(XX, YY, T_norm(:,:,2)), shading flat, colorbar,colormap('cool'), hold on
        % caxis([0 0.12]);
        hold on
        %     dU = U - Uprev;
        %     dV = V - Vprev;
        ppp = quiver( XX(1:1:end,1:1:end), YY(1:1:end,1:1:end), U(1:1:end,1:1:end,1), ...
            V(1:1:end,1:1:end,1), 3, 'k' )
        hold on
        xlabel('x/H'), ylabel('y/H')
        title(['time =', num2str(t)])
        
        %     title(['Thermal contour, t = ', num2str(t)])
        %     xlabel('Top Surface Domain  in the x direction');
        %     ylabel('Top Surface Domain in the y direction')
        pause
    end
   
    
    
    figname = 'Tcontours_fig_';
    if saveFigure == true
        set(ppp, 'PaperPositionMode', 'manual');
        set(ppp, 'PaperUnits', 'inches');
        set(ppp, 'PaperPosition', [0 0 10 9]);
        %        fullfigname = sprintf('%s%d%s',figname, n-1+1000, '.eps');
        fullfigname = sprintf('%s%d%s',figname, t-1+1000, '.jpg');
        %        print(p, '-depsc', fullfigname);
        print(ppp, '-djpeg', fullfigname);
    end
end
