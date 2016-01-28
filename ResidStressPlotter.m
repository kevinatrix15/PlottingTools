%% This script plots calculated residual stress values
clear all, close all, clc

tot_steps = 900;
data_dir = pwd;
% data_dir = 'D:\3DSIM_Source\MyProjects\3DSIMSolver\SolverData\Output_061015';
plot_tool_dir = 'C:\Users\Kevin\OneDrive\3DSIM';
cd(data_dir);

Nodes = 32;
zLayers = 5;

filename = ['Result_T0_ResidStressVals.csv'];
M = csvread(filename, 1);
dx = M(2,1) - M(1,1);
dy = M(Nodes+1,2) - M(1,2);

Lx = Nodes*dx;  Ly = Nodes*dy;

XX = 0.5*dx:dx:(Lx-0.5*dx);
YY = 0.5*dy:dy:(Ly-0.5*dy);

for t=1:1:tot_steps
    filename = ['Result_T', num2str(t-1), '_ResidStressVals.csv'];
    M = csvread(filename, 1);
    
    dz = M(Nodes*Nodes+1, 3) - M(1,3);      % dz varies with time
    Lz = dz*(zLayers-1);
    ZZ = -0.5*Lz:dz:0.5*Lz;
    
    sigxx = zeros(Nodes,Nodes,zLayers);
    sigyy = zeros(Nodes,Nodes,zLayers);
    sigxy = zeros(Nodes,Nodes,zLayers);
    
    count = 1;
    for k = 1:zLayers
        for i = 1:Nodes
            for j = 1:Nodes
                sigxx(i,j,k) = M(count,4);
                sigyy(i,j,k) = M(count,5);
                sigxy(i,j,k) = M(count,6);
                
                count = count+1;
            end
        end
    end
    
%     p = figure(1)
%     %for planeX= 0:1*dz:Lz
%     for zPlane = -0.5*Lz:dz:0.5*Lz;
%         slice(YY, XX, ZZ, sigxx(:,:,:),[],[],zPlane), shading flat, colorbar
%         xlabel('Y'), ylabel('X'), zlabel('Z')
%         hold on
%         pause
%     end
%     %     caxis([297 300]);
%     hold off
    
%     figure(2)
%     for i=1:zLayers
%         subplot(1,5,i)
%         contourf(XX, YY, sigxx(:,:,i))
% %         caxis([min(min(min(sigxx))), max(max(max(sigxx)))])
%         caxis([-1.2e8, 1.2e8])
%         caxis manual;
%         title(['z/t = ', num2str(ZZ(i)/Lz)]);
%         xlabel('X'), ylabel('Y');
% %         if (i==5)
% %             colorbar
% %         end
%         
% %         subplot(2,5,i+5)
% %         contourf(XX, YY, sigxy(:,:,i))
% % %         caxis([-1.2e8, 1.2e8])
% %         caxis manual
% %         title(['z/t = ', num2str(ZZ(i)/Lz)]);
% %         xlabel('X'), ylabel('Y');
% %         if (i==5)
% %             colorbar
% %         end
%         
%     end
%     pause(0.1)
%     
    [XXi, YYi] = meshgrid(XX, YY);
    [XXf, YYf] = meshgrid(0.5*dx:dx/20:(Lx-0.5*dx), 0.5*dy:dy/20:(Ly-0.5*dy));
    
    figure(3)
    subplot(1,2,1)
%     contourf(XX, YY, sigxx(:,:,1))
    contourf(XXf, YYf, interp2(XXi, YYi, squeeze(sigxx(:,:,1)), XXf, YYf))
    h = colorbar
    ylabel(h, '$$\sigma_{xx}, [Pa]$$', 'Interpreter', 'latex', 'FontSize', 14)
    caxis([-3e9, 4e9])
    caxis manual;
    title(['$$z/t = ', num2str(ZZ(1)/Lz), '$$'],'Interpreter','latex', 'FontSize', 14);
    xlabel('$$x$$','Interpreter','latex', 'FontSize', 14), ylabel('$$y$$','Interpreter','latex', 'FontSize', 14);
    
    subplot(1,2,2)


    contourf(XXf, YYf, interp2(XXi, YYi, squeeze(sigxx(:,:,5)), XXf, YYf))
%     contourf(XX, YY, -sigxx(:,:,1))
    h = colorbar
    ylabel(h, '$$\sigma_{xx}, [Pa]$$', 'Interpreter', 'latex', 'FontSize', 14)
    caxis([-3e9, 4e9])
    caxis manual;
    title(['$$z/t = ', num2str(ZZ(5)/Lz), '$$'],'Interpreter','latex', 'FontSize', 14);
    xlabel('$$x$$','Interpreter','latex', 'FontSize', 14), ylabel('$$y$$','Interpreter','latex', 'FontSize', 14);
    
    hold on
    plot([0 2e-3], [2e-3 0], '--k', 'LineWidth', 2), hold on
    plot([0 2e-3], [0 2e-3], '--k', 'LineWidth', 2)
    text(0.2e-3, 1.8e-3, num2str(t))
    hold off
    pause(0.1)
    % 41 (beginning) 455 (end)
    
%     %----- Deformation -----
%     filename2 = ['Result_T', num2str(t-1), '_LocalDisplacements.csv'];
%     N = csvread(filename2, 1);
%     
% %     dz = M(Nodes*Nodes+1, 3) - M(1,3);      % dz varies with time
% %     Lz = dz*(zLayers-1);
% %     ZZ = -0.5*Lz:dz:0.5*Lz;
%     
%     u = zeros(Nodes,Nodes,zLayers);
%     v = zeros(Nodes,Nodes,zLayers);
%     w = zeros(Nodes,Nodes,zLayers);
%     
%     count = 1;
%     for k = 1:zLayers
%         for i = 1:Nodes
%             for j = 1:Nodes
%                 u(i,j,k) = N(count,4);
%                 v(i,j,k) = N(count,5);
%                 w(i,j,k) = N(count,6);
%                 
%                 count = count+1;
%             end
%         end
%     end
%     
%     figure(4)
% %     for i=1:zLayers
%         subplot(1,2,1)
%         contourf(XX, YY, w(:,:,5))
%         colorbar
% %         caxis([min(min(min(u))), max(max(max(u)))])
%         caxis([-6e-6, 0])
%         caxis manual;
%         title(['z/t = ', num2str(ZZ(5)/Lz)]);
%         xlabel('X'), ylabel('Y');
% %         if (i==5)
% %             colorbar
% %         end
%         
%         subplot(1,2,2)
%         contourf(XX, YY, v(:,:,5))
%         caxis([min(min(min(v))), max(max(max(v)))])
% %         caxis([-1.2e8, 1.2e8])
%         caxis manual
%         title(['z/t = ', num2str(ZZ(5)/Lz)]);
%         xlabel('X'), ylabel('Y');
% %         if (i==5)
%             colorbar
% %         end
%         
% %     end
%     pause(0.1)
pause
    
end

                
