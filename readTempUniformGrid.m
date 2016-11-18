function [xvec, yvec, Temps] = ...
  readTempUniformGrid(curr_dir, data_path, filename)

  cd(data_path);

  N = csvread(filename, 1);

  xlocs = N(:,1);
  ylocs = N(:,2);
  zlocs = N(:,3);
  Temps = N(:,4);

  num_pts = length(xlocs);

  nx = find(ylocs(2:end) ~= ylocs(1), 1, 'first');
  ny = num_pts / nx;

  dx = xlocs(2) - xlocs(1);
  dy = ylocs(nx+1) - ylocs(1);

  Lx = xlocs(end);
  Ly = ylocs(end);

  xvec = 0:dx:Lx;
  yvec = 0:dy:Ly;

  cd(curr_dir);
