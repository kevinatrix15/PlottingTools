clear all, clc

plot_tool_dir = '/home/kevin/3DSIM_src/plotting-tools';
saveFigure = false;

%----- Temperature Params ----
plotTemperature = true;
temperatureFileStem = 'frame';
%temperatureFileStem = 'Temperature';
tempFileRoot = '';
% Tfield_dir = '../Output_102715_Animation';
% Tfield_dir= '/home/kevin/3DSIM_src/thermal-solver/SolverData/3dsim-part-output';
% Tfield_dir= '/home/kevin/3DSIM_src/thermal-solver/SolverData/output-1x1x1-fine';
% Tfield_dir= '/home/kevin/3DSIM_src/thermal-solver/SolverData/Output';
% Tfield_dir= '/home/kevin/3DSIM_src/thermal-solver/SolverData/output';
% Tfield_dir = '/home/kevin/3DSIM_src/thermal-solver/test/output';
% Tfield_dir = '/home/kevin/3DSIM_src/thermal-solver/SolverData/gus-0p125x-100um-freq30-layer275';
% Tfield_dir = pwd;
% Tfield_dir = '/home/kevin/Downloads/5x5x5-newTS/10-10/output';
 Tfield_dir = '/home/kevin/Downloads/domframe';
% Tfield_dir = '/home/kevin/sim-output/formnext-vis/lshape/output';
% Tfield_dir = '/home/kevin/Downloads/2x2x2-tuningReheating/25-25/output';
% Tfield_dir = '/home/kevin/Downloads/2x2x2-oldTS/10-10/output';

%----- Strain Parameters -----
% strain_dir = '/home/kevin/3DSIM_src/thermal-solver/SolverData/gus-0p125x-100um-freq30-layer275/shrinkage';
% strain_dir = '/home/kevin/3DSIM_src/thermal-solver/SolverData/gus-0p125x-100um-freq30/shrinkage';
% strain_dir = '/home/kevin/3DSIM_src/thermal-solver/SolverData/xshape-output/shrinkage';
% strain_dir = '/home/kevin/3DSIM_src/thermal-solver/SolverData/output-1x1x1-fine/shrinkage';
% strain_dir = '/home/kevin/3DSIM_src/thermal-solver/SolverData/output/shrinkage';
% strain_dir = '/home/kevin/3DSIM_src/thermal-solver/SolverData/Output/shrinkage';
% strain_dir = '/home/kevin/3DSIM_src/thermal-solver/test/reference-linux-periodic/shrinkage';
% strain_dir = '/home/kevin/3DSIM_src/thermal-solver/test/output-periodic/shrinkage';
% strain_dir = '/home/kevin/3DSIM_src/thermal-solver/SolverData/3dsim-part-output/shrinkage';
% strain_dir = '/home/kevin/3DSIM_src/materialise-prep/thermalsimulation-319/output/shrinkage';
% strain_dir = '/home/kevin/3DSIM_src/materialise-prep/thermalsimulation-269/output/shrinkage';
% strain_dir = '/home/kevin/Downloads/shrinkage_v293'; % Failed
% strain_dir = '/home/kevin/Downloads/shrinkage_253_oshape';
% strain_dir = '/home/kevin/Downloads/shrinkage_260';
% strain_dir = '/home/kevin/Downloads/shrinkage_298_longcool';
% strain_dir = '/home/kevin/Downloads/shell-0p125x-shrinkage';
% strain_dir = '/home/kevin/Downloads/gus-0p125x-200um-freq50';
% strain_dir = '/home/kevin/Downloads/shrinkage';
% strain_dir = '/home/kevin/Downloads/shrinkage_297_shortcool';
% strain_dir = '/home/kevin/Downloads/shrinkage_264';
strain_dir = '/home/kevin/Downloads/5x5x5-newTS/100-100/output/shrinkage';
% strain_dir = '/home/kevin/Downloads/2x2x2-oldTS/10-10/output/shrinkage';
% strain_dir = '/home/kevin/Downloads/shell-0p125x-50um-20freq-2p26';
% strain_dir = [Tfield_dir, '/shrinkage'];

plotStrain = false;
shrinkageFileStem = 'Shrinkage';
% shrinkageFileStem = 'ShrinkageDisplacementsShrinkage';


%----- Mesh Parameters -----
sd=12;   % number of square divisions in the mesh
Nodes=21;

cd(plot_tool_dir)
 [X, Y, ~] = readTempUniformGrid(pwd, Tfield_dir, [tempFileRoot, ...
     temperatureFileStem, '.csv']);
% [X, Y, ~] = readStrain(pwd, strain_dir, [shrinkageFileStem, '_L10001.csv']);
cd(strain_dir)

xlocs_full = zeros(length(X)*length(Y), 1);
ylocs_full = zeros(length(X)*length(Y), 1);
for j = 1:length(Y)
  for i = 1:length(X)
    idx = i + (j-1)*length(X);
    xlocs_full(idx) = X(i);
    ylocs_full(idx) = Y(j);
  end
end

max_x= X(end);
min_x=0;
max_y= Y(end);
min_y=0;

bedwidth_x=1.2*(max_x-min_x);
bedwidth_y=1.2*(max_y-min_y);
bedmin_x=min_x;
bedmax_x=max_x+0.2*bedwidth_x;
bedmin_y=min_y;
bedmax_y=max_y+0.2*bedwidth_y;


%----- Set up grid to interpolate to -----
[qx,qy] = meshgrid([bedmin_x:bedwidth_x/(sd*(Nodes-1)):bedmax_x],[bedmin_y:bedwidth_y/(sd*(Nodes-1)):bedmax_y]);

%----- Read all strain file names and sort -----
if (plotStrain == true)
    cd(strain_dir)
    strainFiles = dir([shrinkageFileStem, '*']);
    [~,index] = sort_nat({strainFiles.name});
    strainFiles = strainFiles(index);
    numFiles = length(strainFiles);
end

%----- Read the first layer of temperature -----
if (plotTemperature == true)
    cd(Tfield_dir)

    Tfiles = dir(['*', temperatureFileStem, '*']);
    [~,index] = sort_nat({Tfiles.name});
    Tfiles = Tfiles(index);
    numFiles = length(Tfiles);

    % set(gca, 'nextplot', 'replacechildren');
    [qx,qy] = meshgrid([bedmin_x:bedwidth_x/(sd*(Nodes-1)):bedmax_x],[bedmin_y:bedwidth_y/(sd*(Nodes-1)):bedmax_y]);
%     [~,index] = sortrows([Tfiles.datenum].'); Tfiles = Tfiles(index);
    cd(strain_dir)
end

%----- Loop over all files ----
for t = 1:1:numFiles
    clf
    ppp = figure(10)

    if (plotTemperature == true)
        cd(Tfield_dir)
        SallyT=csvread(Tfiles(t).name);
        cd(strain_dir)
        p=SallyT(:,1:2);
        %     plot(p(:,1),p(:,2),'r*')
        T2=SallyT(:,4);
        F = TriScatteredInterp(p,T2);
        qz = F(qx,qy);
        %             mesh(qx,qy,qz);


        hold on
         % subplot(1,2,1)
        contourf(qx,qy,qz, [350 400 500 600 700 1000 1200 1500 1800 1923 3315]);
        caxis([350, 3315]);
        caxis manual;
        %              caxis([(300)/max(max(T)) (1500)/max(max(T))]*max(max(T)));

        colorbar
        title(['time =', num2str(t)])
        % axis([0 0.0215 0 0.0115])
        % axis([0.65e-3 2.6e-3 0.65e-3 2.6e-3])
        axis([0 0.041 0 0.011])
    end

    %----- Plot Strain -----
    if plotStrain == true
        filename = strainFiles(t).name;
        % filename = [shrinkageFileStem, '_L', num2str(t+10000), '.csv'];

        cd(plot_tool_dir)
        [X, Y, Eps] = readStrain(pwd, strain_dir, filename);
        cd(strain_dir)

        F = TriScatteredInterp(xlocs_full, ylocs_full, Eps);
        qz = F(qx,qy);
         % subplot(1,2,2)
        ppp = contourf(qx,qy,qz);

        % axis([0 0.011 0 0.003])
        % axis([0 0.0215 0 0.0115])

        hold on
        xlabel('x'), ylabel('y')
        %title(['layer =', num2str(t)])
        title(['time =', num2str(t)])

        caxis([-0.02, 0.02])
        colorbar
        colormap(jet)
        pause(0.005)
        hold on
    end

    pause(0.001)


    figname = 'strain_';
    % figname = 'strain_temp_';
    if saveFigure == true
        % set(ppp, 'PaperPositionMode', 'manual');
        % set(ppp, 'PaperUnits', 'inches');
        % set(ppp, 'PaperPosition', [0 0 10 5]);
        %        fullfigname = sprintf('%s%d%s',figname, n-1+1000, '.eps');
        fullfigname = sprintf('%s%d%s',figname, t-1+10000, '.jpg');
        %        print(p, '-depsc', fullfigname);
        % print(ppp, '-djpeg', fullfigname);
        saveas(gca, fullfigname);
    end
end
