function rot = getHeadingRotatedVector(heading, vec)
% heading   : Nx1 heading (rad)
% vec       : 1x3 mag vector (micro tesla)
% rot       : Nx3 mag vector (micro tesla)

% Create rotate matrix depending on magnetic heading, in each particles
% rot = arrayfun(@(x)(([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]*vec')'),heading,'UniformOutput',false);
% rot = cell2mat(rot);
% rot = rot';

mag_x = arrayfun(@(x) cos(x)*vec(1)-sin(x)*vec(2), heading);
mag_y = arrayfun(@(x) sin(x)*vec(1)+cos(x)*vec(2), heading);
rot = [mag_x, mag_y, vec(3)*ones(size(heading))];