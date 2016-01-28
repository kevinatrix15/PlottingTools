%% Read in Surface Roughness .csv file
clear all, clc

nx = 25;
ny = 25;

npx = 85*nx;
npy = 85*ny;
dx = 0.0025/npx;
dy = 0.0025/npy;


powers = [310 335 360 385];
speeds = [850 1000 1200 1400 1600 1800 2000 2250 2500];%

% Setup output matrix
outMat = zeros(length(powers)*length(speeds)+1, 5);
% outMat(1,1) = 'Power';
% outMat(1,2) = 'Speed';
% outMat(1,3) = 'Mean Height (m)';
% outMat(1,4) = 'Variance';
% outMat(1,5) = 'STD (m)';



row = 1;
for p = 1:length(powers)
    for s = 1:length(speeds)
        filename = ['Result_TopSurfaceRoughness.csv']
        %         filename = ['GE_SolverResult_', num2str(powers(p)),...
        %             '_', num2str(speeds(s)), '.csv_TopSurfaceRoughness.csv']
        M = csvread(filename, 4, 0);
        xp = M(:,1);
        yp = M(:,2);
        zp = M(:,3);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assign values from file
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j=1:npy
            for i=1:npx
                X(i,j) = dx*(i+0.5);
                Y(i,j) = dy*(j+0.5);
                Z(i,j) = zp(i + (j-1)*npx);
            end
        end
       
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute surface height statistics
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Mean Height
        count = 0;
        sum = 0;
        for j=1:npy
            for i=1:npx
                if ((Y(i,j) < 0.0021 && Y(i,j) > 0.00035) && ...
                        (X(i,j) > 0.00025 && X(i,j) < 0.00225))
                    sum = sum + Z(i,j);
                    count = count + 1;
                end
            end
        end
        avg_height = sum/count
        
        % Variance
        var_sum = 0;
        for j=1:npy
            for i=1:npx
                if ((Y(i,j) < 0.0021 && Y(i,j) > 0.00035) && ...
                        (X(i,j) > 0.00025 && X(i,j) < 0.00225))
                    var_sum = var_sum + (Z(i,j)-avg_height)^2;
                end
            end
        end
        variance = var_sum/count
        std = sqrt(variance)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Store values to matrix
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        row = row+1;
        outMat(row,1) = powers(p);
        outMat(row,2) = speeds(s);
        outMat(row,3) = avg_height;
        outMat(row,4) = variance;
        outMat(row,5) = std;

       
    end
end


for p = 1:length(powers)
    begin_row = (p-1)*length(speeds) + 2;
    end_row = begin_row+length(speeds) - 1;
    
    % Mean height
    figure(1)
    plot(outMat(begin_row:end_row,2), outMat(begin_row:end_row,3))
    hold on
    
    % STD
    figure(2)
    plot(outMat(begin_row:end_row,2), outMat(begin_row:end_row,5))
    hold on
end




%% Read in Surface Roughness .csv file
clear all, clc

% M = csvread('GE_SolverResult_310_1800.csv_TopSurfaceRoughness.csv', 4, 0);
% M = csvread('Result_testbed_TopSurfaceRoughness.csv', 4, 0);
M = csvread('Result_TopSurfaceRoughness.csv', 4, 0);


nx = 25;
ny = 25;

xp = M(:,1);
yp = M(:,2);
zp = M(:,3);
npx = 85*nx;
npy = 85*ny;
dx = 0.0025/npx;
dy = 0.0025/npy;

for j=1:npy
    for i=1:npx
        X(i,j) = dx*(i+0.5);
        Y(i,j) = dy*(j+0.5);
        Z(i,j) = zp(i + (j-1)*npx);
    end
end

h=figure(2)
surf(X,Y,Z, 'EdgeColor', 'none','LineStyle','none','FaceLighting','gouraud')
axis([0 0.0025 0 0.0025 0 2e-4])
colorbar, caxis([0 0.5e-4])
set(h, 'Color', [1 1 1])


% Compute surface height statistics

% Mean
count = 0;
sum = 0;
for j=1:npy
    for i=1:npx
        if ((Y(i,j) < 0.0021 && Y(i,j) > 0.00035) && ...
                (X(i,j) > 0.00025 && X(i,j) < 0.00225))
            sum = sum + Z(i,j);
            count = count + 1;
        end
    end
end

avg_height = sum/count

% Variance
var_sum = 0;
for j=1:npy
    for i=1:npx
        if ((Y(i,j) < 0.0021 && Y(i,j) > 0.00035) && ...
                (X(i,j) > 0.00025 && X(i,j) < 0.00225))
            var_sum = var_sum + (Z(i,j)-avg_height)^2;
        end
    end
end
variance = var_sum/count

%% 
clear all, clc
N = csvread('Result_L3_Meltpool.csv', 1, 0);

tstep = N(:,3);
midX = N(:,18);
midY = N(:,19);

% xx = linspace(min(midX):
xx = linspace(0,0.003,10)
yy = linspace(0,0.003,10)

figure(2)
plot(midX, midY, 'o')
axis([0 max(xx) 0 max(yy)])



%% 
[X,Y,Z] = peaks(25)
figure
surf(X,Y,Z)


