function [U, V, Eps, xvec, yvec] = readShrinkage(curr_dir, data_path, filename)

  cd(data_path);

  N = csvread(filename, 1);

  xlocs = N(:,1);
  ylocs = N(:,2);
  xdisp = N(:,3);
  ydisp = N(:,4);
  strain = N(:,5);

  num_pts = length(xlocs);

  nx = find(ylocs(2:end) ~= ylocs(1), 1, 'first');
  ny = num_pts / nx;

  dx = xlocs(2) - xlocs(1);
  dy = ylocs(nx+1) - ylocs(1);

  Lx = xlocs(end);
  Ly = ylocs(end);

  xvec = 0:dx:Lx;
  yvec = 0:dy:Ly;

  U = zeros(nx,ny);
  V = zeros(nx,ny);
  Eps = zeros(nx,ny);
  for i= 1:nx
      for j=1:ny
          U(i,j) = xdisp(i + (j-1)*nx);
          V(i,j) = ydisp(i + (j-1)*nx);
          Eps(i,j) = strain(i + (j-1)*nx);
      end
  end

  cd(curr_dir);
