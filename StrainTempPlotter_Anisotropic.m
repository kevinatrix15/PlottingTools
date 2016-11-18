clear all, clc

plot_tool_dir = '/home/kevin/3DSIM_src/plotting-tools';
saveFigure = false;

%----- Temperature Params ----
plotTemperature = false;
temperatureFileStem = 'Temperature';
 Tfield_dir = '/home/kevin/Downloads/2x2x2-tuningReheating/25-25/output';
% Tfield_dir = '/home/kevin/3DSIM_src/thermal-solver/test/output';

%----- Strain Parameters -----
% strain_dir = '/home/kevin/3DSIM_src/thermal-solver/SolverData/output/shrinkage';
% strain_dir = '/home/kevin/3DSIM_src/thermal-solver/test/reference-linux-periodic/shrinkage';
% strain_dir = '/home/kevin/3DSIM_src/thermal-solver/test/output/shrinkage';
% strain_dir = '/home/kevin/Downloads/5x5x5-beta/10-10/output/shrinkage';
% strain_dir = '/home/kevin/Downloads/5x5x5-threshold/25-25/output/shrinkage';
% strain_dir = '/home/kevin/Downloads/2x2x2-newTS/25-25/output/shrinkage';
% strain_dir = '/home/kevin/Downloads/2x2x2-tuningReheating/25-25/output/shrinkage';
% strain_dir = '/home/kevin/sim-output/debug/vertplate/strain';
 strain_dir = '/home/kevin/sim-output/debug/aibuild-1661/output/shrinkage';
% strain_dir = '/home/kevin/Downloads/2x2x2-4.16.0/100-100/output/shrinkage';
%strain_dir = '/home/kevin/3DSIM_src/thermal-solver/test/reference-linux-periodic/shrinkage';
% strain_dir = '/home/kevin/3DSIM_src/thermal-solver/test/output_20e6_mf1_vtk/shrinkage';
% strain_dir = '/home/kevin/3DSIM_src/thermal-solver/test/output/shrinkage';
% strain_dir = '/home/kevin/Downloads/shrinkage';
% strain_dir = '/home/kevin/Downloads/2x2x2/25-25/output/shrinkage';
%strain_dir = [Tfield_dir + '/shrinkage'];

plotStrain = true;
shrinkageFileStem = 'Shrinkage';

YS = 5.4e8;
E = 1.96e11;
rf = 2e-2;
yieldStrain = 1;
% yieldStrain = rf*YS/E;

%----- Mesh Parameters -----
sd=12;   % number of square divisions in the mesh
Nodes=21;

cd(plot_tool_dir)
  if (plotTemperature == true && plotStrain == false)
    [X, Y, ~] = readTempUniformGrid(pwd, Tfield_dir, ['ResultL0_T0_', ...
        temperatureFileStem, '.csv']);
  else
    [X, Y, ~, ~, ~] = readAnisotropicStrain(pwd, strain_dir, [shrinkageFileStem, '_L10001.csv']);
  end
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

    cd(plot_tool_dir)
    [~,index] = sort_nat({strainFiles.name});
    strainFiles = strainFiles(index);
    numFiles = length(strainFiles);
    cd(strain_dir)
end

%----- Read the first layer of temperature -----
if (plotTemperature == true)
    cd(Tfield_dir)

    Tfiles = dir(['*', temperatureFileStem, '*']);
    cd(plot_tool_dir)
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
    ppp = figure(10)

    if (plotTemperature == true)
        cd(Tfield_dir)
        SallyT=csvread(Tfiles(t).name);

        cd(plot_tool_dir)
        [X, Y, ~] = readTempUniformGrid(pwd, Tfield_dir, Tfiles(t).name);
        cd(strain_dir)

        p=SallyT(:,1:2);
        %     plot(p(:,1),p(:,2),'r*')
        T2=SallyT(:,4);
        F = TriScatteredInterp(p,T2);
        qz = F(qx,qy);
        %             mesh(qx,qy,qz);


        hold on
        contourf(qx,qy,qz, [350 400 500 600 700 1000 1200 1500 1800 1923 3315]);
        caxis([350, 3315]);
        caxis manual;
        %              caxis([(300)/max(max(T)) (1500)/max(max(T))]*max(max(T)));

        colorbar
        title(Tfiles(t).name);
        %axis([min(X) 0.002 min(Y) 0.002])
        % axis([0 0.0215 0 0.0115])
    end

    %----- Plot Anisotropic Strain -----
    if plotStrain == true
        filename = strainFiles(t).name;
        % filename = [shrinkageFileStem, '_L', num2str(t+10000), '.csv'];

        cd(plot_tool_dir)
        [X, Y, EpsXX, EpsYY, EpsZZ, EpsXY] = readAnisotropicStrain(pwd, strain_dir, filename);
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

        Fxx = TriScatteredInterp(xlocs_full, ylocs_full, EpsXX);
        qz_xx = Fxx(qx,qy);
         subplot(2,2,1)
        ppp = contourf(qx,qy,qz_xx/yieldStrain);

        title('\epsilon_x_x', 'fontsize', 14)
        %title(['time =', num2str(t)])
        % axis([0 1.2e-3 0 1.2e-3])
        % caxis([-1e-3, 0])
        % caxis([-3, 0])
        colorbar
        colormap(jet)

        Fyy = TriScatteredInterp(xlocs_full, ylocs_full, EpsYY);
        qz_yy = Fyy(qx,qy);
         subplot(2,2,2)
        ppp = contourf(qx,qy,qz_yy/yieldStrain);

        title('\epsilon_y_y', 'fontsize', 14)
        % axis([0 1.2e-3 0 1.2e-3])
        %caxis([-0.02, 0.007])
         %caxis([-1e-3, 0])
        % caxis([-3, 0])
        colorbar
        colormap(jet)

        Fzz = TriScatteredInterp(xlocs_full, ylocs_full, EpsZZ);
        qz_zz = Fzz(qx,qy);
        subplot(2,2,3)
        ppp = contourf(qx,qy,qz_zz/yieldStrain);

        title('\epsilon_z_z', 'fontsize', 14)
        % axis([0 1.2e-3 0 1.2e-3])
        % caxis([-0.02, 0.007])
         %caxis([-3e-4, 0])
        % caxis([-1, 0])
        colorbar
        colormap(jet)

        Fxy = TriScatteredInterp(xlocs_full, ylocs_full, EpsXY);
        qz_xy = Fxy(qx,qy);
         subplot(2,2,4)
        ppp = contourf(qx,qy,qz_xy/yieldStrain);

        % axis([0 0.011 0 0.003])
        % axis([0 0.0215 0 0.0115])

        hold on
        xlabel('x'), ylabel('y')
        %title(['layer =', num2str(t)])
        title(['time =', num2str(t)])

        % title('\epsilon_x_y', 'fontsize', 14)
        % axis([0 1.2e-3 0 1.2e-3])
        % caxis([-0.02, 0.007])
         %caxis([-4e-4, 0])
        % caxis([-1, 0])
        colorbar
        colormap(jet)
        pause(0.005)
        hold on
    end

    pause


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
