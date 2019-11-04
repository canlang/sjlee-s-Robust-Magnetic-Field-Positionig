clear; close all;
addpath('xyz file operations')

%%
data1 = readtable('input_traj/result.csv');


lM = [data1.mag_x,data1.mag_y,data1.mag_z];
% DATA SPLIT
x = data1.gt_x;
y = data1.gt_y;
% z = vecnorm(lM,2,2);    % L2 norm
z = data1.mag_x;     % x component
% z = data1.mag_y;     % y component
% z = data1.mag_z;     % z component

tranR = arrayfun(@(a,b,c,d,e,f) euler2rotMat(a,b,c)*[d;e;f]...
    ,data1.roll,data1.pitch,data1.yaw,lM(:,1),lM(:,2),lM(:,3),'un',0);
tranR2 = reshape(cell2mat(tranR),3,[]).';
z = tranR2(:,1);
lM = tranR2;
data1.x = data1.gt_x;
data1.y = data1.gt_y;
%%
A = imread('N1-7F.png','BackgroundColor',[1 1 1]);
xWorldLimits = [-1 1650/20];
yWorldLimits = [-1 660/20];
RA = imref2d(size(A),xWorldLimits,yWorldLimits);
% imshow(flipud(A),RA);

hold on
[X,Y,Z] = interpolation_by_alphashape(x,y,z,0.1);

h = contourf(X,Y,Z);
hcb = colorbar;
hold off
% %% ----------------------------------------------------
% xlabel('(m)')
% ylabel('y (m)')
% axis image
% axis xy
% 
% set(gcf,'units','points','position',[700,500,1000,350])
% tightfig
% legend('alphashape')
set(gcf,'units','points','position',[500,500,1600,600])
sdf(gcf,'sj2')


%% INTERPOLATION by alphashape
function [X,Y,Z] = interpolation_by_alphashape(x,y,z,interval)
XI = min(x):interval:max(x);
YI = min(y):interval:max(y);
[X,Y] = meshgrid(XI,YI);

% APPLY ALPHA SHAPE 
shp = alphaShape(x,y);
% plot(shp,'FaceAlpha',0.5,'EdgeAlpha',0.3,'LineWidth',0.3)
plot(shp,'FaceColor','red','FaceAlpha',0.3,'LineWidth',0.3)
% 'FaceColor','red',

% hline = findobj(gcf, 'type', 'line');
% set(hline(1),'LineStyle','--')
% shp.Alpha = 2;            % Alpha parameter: 
in = inShape(shp,X,Y);
xg = X(in);
yg = Y(in);
% DATA INPUT ONLY INSIDE OF ALPHA SHAPE
zg = griddata(x,y,z,xg,yg);         % 2. griddata() : INTERPOLATION
[X,Y,Z] = xyz2grid(xg,yg,zg);
end
