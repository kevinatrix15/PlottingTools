theta = -45;
qpx = sqrt(2)/2;
qpy = 0;
p0x = 0;
p0y = 0;

% qx = (qpx + p0x)*cosd(theta) - (qpy + p0y)*sind(theta)
% qy = (qpx + p0x)*sind(theta) + (qpy + p0y)*cosd(theta)

qx = (qpx)*cosd(theta)+ p0x - (qpy)*sind(theta)
qy = (qpx)*sind(theta) + (qpy)*cosd(theta) + p0y

%%
theta = 30;
qx = sqrt(3)/2;
qy = 0.5;
p0x = 0;
p0y = 0;

qpx = (qx - p0x)*cosd(theta) + (qy - p0y)*sind(theta)
qpy = (p0x - qx)*sind(theta) + (qy - p0y)*cosd(theta)
% qx = (qpx + p0x)*cosd(theta) - (qpy + p0y)*sind(theta)
% qy = (qpx + p0x)*sind(theta) + (qpy + p0y)*cosd(theta)
