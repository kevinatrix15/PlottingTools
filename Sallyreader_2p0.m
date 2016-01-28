N=1000;
windowing=1;
windowing1=20; 
Nodes=21;
sd=4;

temperatureFileStem = 'thermal';

mpfile=csvread([temperatureFileStem, '_L1_Meltpool.csv'],1,0);
% max_x=max(mpfile(:,4));
% min_x=min(mpfile(:,4));
% max_y=max(mpfile(:,5));
% min_y=min(mpfile(:,5));
max_x=16e-3;
min_x=0;
max_y=6e-3;
min_y=0;

bedwidth_x=1.2*(max_x-min_x);
bedwidth_y=1.2*(max_y-min_y);
bedmin_x=min_x;
bedmax_x=max_x+0.2*bedwidth_x;
bedmin_y=min_y;
bedmax_y=max_y+0.2*bedwidth_y;

SallyT=csvread([temperatureFileStem, '_L1_O1_T0_Temperature.csv']);

% SallyT=csvread('Result_L1_O1_T0_Temperature.csv');
% set(gca, 'nextplot', 'replacechildren');
[qx,qy] = meshgrid([bedmin_x:bedwidth_x/(sd*(Nodes-1)):bedmax_x],[bedmin_y:bedwidth_y/(sd*(Nodes-1)):bedmax_y]);
csvfiles = dir('*Temperature.csv'); 
% filenames = cellstr(csvfiles(:).name);
% [filenames, idx] = sort_nat(filenames, 'ascend');
[~,index] = sort_nat({csvfiles.name}); 
csvfiles = csvfiles(index);

% [~,index] = sortrows([csvfiles.datenum].'); csvfiles = csvfiles(index);
numfiles = length(csvfiles);

saveFigure = false;

for i = 3650:1:numfiles
            clf
            ppp = figure(1)
            SallyT=csvread(csvfiles(i).name);
            p=SallyT(:,1:2);
%             plot(p(:,1),p(:,2),'k*')
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
            title(['Thermal contour, t = ', num2str(i)])
            xlabel('Top Surface Domain  in the x direction');
            ylabel('Top Surface Domain in the y direction')
%             pause(0.001)
            pause

            
            figname = 'Tcontours_fig_';
            if saveFigure == true
                set(ppp, 'PaperPositionMode', 'manual');
                set(ppp, 'PaperUnits', 'inches');
                set(ppp, 'PaperPosition', [0 0 13 6]);
                %        fullfigname = sprintf('%s%d%s',figname, n-1+1000, '.eps');
                fullfigname = sprintf('%s%d%s',figname, i-1+1000, '.jpg');
                %        print(p, '-depsc', fullfigname);
                print(ppp, '-djpeg', fullfigname);
            end
end
