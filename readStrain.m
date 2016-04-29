function [xvec, yvec, strain] = readShrinkage(curr_dir, data_path, filename)

  cd(data_path);

  N = csvread(filename, 1);

  xlocs = N(:,1);
  ylocs = N(:,2);
  strain = N(:,3);

  num_pts = length(xlocs);

  nx = find(ylocs(2:end) ~= ylocs(1), 1, 'first');
  ny = num_pts / nx;

  dx = xlocs(2) - xlocs(1);
  dy = ylocs(nx+1) - ylocs(1);

  Lx = xlocs(end);
  Ly = ylocs(end);

  xvec = 0:dx:Lx;
  yvec = 0:dy:Ly;

%  Eps = zeros(nx,ny);
%  for i= 1:nx
%      for j=1:ny
%          Eps(i,j) = strain(i + (j-1)*nx);
%      end
%  end

  cd(curr_dir);
