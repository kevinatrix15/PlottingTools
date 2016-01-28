function [x, y, vel] = readKeyholeVels(filename, data_dir)

curr_dir = pwd;
cd(data_dir)
if exist(filename, 'file')
    M = csvread(filename, 1);
else
    M = zeros(1,3);
end
cd(curr_dir)

x = M(:,1);
y = M(:,2);
vel = M(:,3);

return;
end

