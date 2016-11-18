function [xvec, yvec, strainXX, strainYY, strainZZ, strainXY] = readAnisotropicStrain(curr_dir, data_path, filename)

  cd(data_path);

  N = csvread(filename, 1);

  xlocs = N(:,1);
  ylocs = N(:,2);
  strainXX = N(:,3);
  strainYY = N(:,4);
  strainZZ = N(:,5);
  strainXY = N(:,6);

  num_pts = length(xlocs);

  nx = find(ylocs(2:end) ~= ylocs(1), 1, 'first');
  ny = num_pts / nx;

  dx = xlocs(2) - xlocs(1);
  dy = ylocs(nx+1) - ylocs(1);

  Lx = xlocs(end);
  Ly = ylocs(end);

  xvec = unique(xlocs);
  yvec = unique(ylocs);
  % xvec = xlocs(1):dx:Lx;
  % yvec = ylocs(1):dy:Ly;

  cd(curr_dir);
