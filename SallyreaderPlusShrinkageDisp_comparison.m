clear all, clc


%----- Set plotting parameters -----
plotTemperature = false;
saveFigure = false;
num_layers = 1;             % Number of layers considered by diagnostic
num_times = 30000;

%----- Shrinkage Parameters -----
% shrinkageFileStem = 'shrinkagefile';
shrinkageFileStem = 'ShrinkageDisplacements';

data_dir1 = 'D:\3DSIM_Source\MyProjects\Solver\SolverData\Output_ss100\shrinkage';
data_dir2 = 'D:\3DSIM_Source\MyProjects\Solver\SolverData\Output_ss100_startLines\shrinkage';
plot_tool_dir = 'C:\Users\Kevin\OneDrive\3DSIM';


Tfield_dir = '../Output_102715_Animation';
% Tfield_dir = pwd;


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
    SallyT=csvread('Result_L1_O1_T0_Temperature.csv');
    
    % set(gca, 'nextplot', 'replacechildren');
    [qx,qy] = meshgrid([bedmin_x:bedwidth_x/(sd*(Nodes-1)):bedmax_x],[bedmin_y:bedwidth_y/(sd*(Nodes-1)):bedmax_y]);
    csvfiles = dir('*Temperature.csv');
    [~,index] = sortrows([csvfiles.datenum].'); csvfiles = csvfiles(index);
    numfiles = length(csvfiles);
    cd(data_dir)
end



for l = 1:num_layers
    
    %----- Temperature Contour Plotting -----
    ppp = figure(1)
    
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
    end
    


    %----- Shrinkage -----
    filename = [shrinkageFileStem, '_L', num2str(l+1000), '.csv'];

    [U1, V1, Eps1, X1, Y1] = readShrinkage(pwd, data_dir1, filename);
    [U2, V2, Eps2, X2, Y2] = readShrinkage(pwd, data_dir2, filename);

    [YY,XX] = meshgrid(Y1, X1);     % should be same size if both grids are equal
    
    subplot(1,2,1)
    quiver( XX(1:1:end,1:1:end), YY(1:1:end,1:1:end), U1(1:1:end,1:1:end,1), ...
        V1(1:1:end,1:1:end,1), 3, 'k' )
    xlabel('x'), ylabel('y')

    subplot(1,2,2)
    quiver( XX(1:1:end,1:1:end), YY(1:1:end,1:1:end), U2(1:1:end,1:1:end,1), ...
        V2(1:1:end,1:1:end,1), 3, 'k' )
    xlabel('x'), ylabel('y')

    %----- Difference between results -----
    errU = (U2 - U1)/max([max(max(U1)) abs(min(min(U1)))]);
    errV = (V2 - V1)/max([max(max(V1)) abs(min(min(V1)))]);
    
    figure(2)
    quiver( XX(1:1:end,1:1:end), YY(1:1:end,1:1:end), U2-U1, V2-V1, 3, 'b' )
    xlabel('x'), ylabel('y')
    
    U1 = reshape(U1,size(U1,1)*size(U1,2),1);
    U2 = reshape(U2,size(U2,1)*size(U2,2),1);
    RMSE = sqrt(sum((U1(:)-U2(:)).^2)/numel(U1))
    
    %----- Save Figures for animation -----
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

% max(max(errU))
% min(min(errU))
