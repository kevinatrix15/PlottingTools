function [x_liq, z_liq, x_vap, z_vap, x_lqcor, z_lqcor, x_vapcor, z_vapcor,...
        x_bkvel, z_bkvel, v_bkvel, x_frvel, z_frvel, v_frvel] = readKeyholeFiles(filename, data_dir)

curr_dir = pwd;
cd(data_dir)
if exist(filename, 'file')
    M = csvread(filename, 2);
else
    M = zeros(1,14);
end
cd(curr_dir)

r_tmp = find(M(:,1)==0,1) - 1;
if isempty(r_tmp)
    r_tmp = size(M, 1)
end
x_liq = M(1:r_tmp,1);
z_liq = M(1:r_tmp,2);

r_tmp = find(M(:,3)==0,1) - 1;
if isempty(r_tmp)
    r_tmp = size(M, 1)
end
x_vap = M(1:r_tmp,3);
z_vap = M(1:r_tmp,4);

r_tmp = find(M(:,5)==0,1) - 1;
if isempty(r_tmp)
    r_tmp = size(M, 1)
end
x_lqcor = M(1:r_tmp,5);
z_lqcor = M(1:r_tmp,6);

r_tmp = find(M(:,7)==0,1) - 1;
if isempty(r_tmp)
    r_tmp = size(M, 1)
end
x_vapcor = M(1:r_tmp,7);
z_vapcor = M(1:r_tmp,8);

r_tmp = find(M(:,9)==0,1) - 1;
if isempty(r_tmp)
    r_tmp = size(M, 1)
end
x_bkvel = M(1:r_tmp,9);
z_bkvel = M(1:r_tmp,10);
v_bkvel = M(1:r_tmp,11);

r_tmp = find(M(:,12)==0,1) - 1;
if isempty(r_tmp)
    r_tmp = size(M, 1)
end
x_frvel = M(1:r_tmp,12);
z_frvel = M(1:r_tmp,13);
v_frvel = M(1:r_tmp,14);

return;
end

