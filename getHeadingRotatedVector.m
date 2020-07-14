function rot = getHeadingRotatedVector(heading, v, R)
% heading   : Nx1 heading (rad)
% v(vec)    : 1x3 mag vector (micro tesla)
% rot       : Nx3 mag vector (micro tesla)

% Create rotate matrix depending on magnetic heading, in each particles
% rot = arrayfun(@(x)(([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]*vec')'),heading,'UniformOutput',false);
% rot = cell2mat(rot);
% rot = rot';

% mag_x = arrayfun(@(x) cos(x)*v(1)-sin(x)*v(2), heading);
% mag_y = arrayfun(@(x) sin(x)*v(1)+cos(x)*v(2), heading);
% rot = [mag_x, mag_y, v(3)*ones(size(heading))];

A = R.';
% A = R;

mag_x = arrayfun(@(x) (cos(x)*A(1)-sin(x)*A(2))*v(1)+(cos(x)*A(4)-sin(x)*A(5))*v(2)+(cos(x)*A(7)-sin(x)*A(8))*v(3), heading);
mag_y = arrayfun(@(x) (sin(x)*A(1)+cos(x)*A(2))*v(1)+(sin(x)*A(4)+cos(x)*A(5))*v(2)+(sin(x)*A(7)+cos(x)*A(8))*v(3), heading);
mag_z = [A(3) A(6) A(9)]*v'*ones(size(heading));
rot = [mag_x,mag_y,mag_z];
% mag_x = arrayfun(@(x) 