clear; close all;

addpath('xyz file operations')

data1 = readtable('batch.csv');

lM = [data1.magnet_x,data1.magnet_y,data1.magnet_z];
% DATA SPLIT
x = data1.x;
y = data1.y;
z = vecnorm(lM,2,2);    % L2 norm
% z = data1.magnet_x;     % x component
% z = data1.magnet_y;     % y component
% z = data1.magnet_z;     % z component


A = imread('N1-7F.png','BackgroundColor',[1 1 1]);
xWorldLimits = [-1 1650/20];
yWorldLimits = [-1 660/20];
RA = imref2d(size(A),xWorldLimits,yWorldLimits);
% imshow(flipud(A),RA);

hold on
[X,Y,Z] = interpolation_by_alphashape(x,y,z);

% h = imagesc(X(1,:),Y(:,1),Z);
% xlabel('x (m)')
% ylabel('y (m)')
% axis image
% axis xy
% set(h,'alphadata',~isnan(Z))
% colormap(ax2,cool)
% colorbar
% F = scatteredInterpolant(data1.magnet_x,data1.magnet_y,data1.magnet_z,'nearest');
% vq = F(X,Y);

%%
% h = contourf(X,Y,Z);
% hcb = colorbar;
% ylabel(hcb, 'L2 norm of Magnetic Field')
% ylabel(hcb, 'X compoment of Magnetic Field')
% ylabel(hcb, 'Y compoment of Magnetic Field')
% ylabel(hcb, 'Z compoment of Magnetic Field')


%% ----------------------------------------------------
xlabel('(m)')
ylabel('y (m)')
axis image
axis xy

set(gcf,'units','points','position',[700,500,1000,350])
tightfig
legend('alphashape')
% set(gcf,'units','points','position',[500,500,1200,600])
sdf(gcf,'sj4')

%% ----------------------------------------------------
% saveas(gcf,'eps/interp_norm.eps','epsc');
% saveas(gcf,'eps/interp_x.eps','epsc');
% saveas(gcf,'eps/interp_y.eps','epsc');
% saveas(gcf,'eps/interp_z.eps','epsc');
% saveas(gcf,'eps/alphashape.eps','epsc');
print -depsc2 eps/alphashape.eps
% export_fig eps/alphashape.eps -depsc -m2.5
%% INTERPOLATION by alphashape
function [X,Y,Z] = interpolation_by_alphashape(x,y,z)
XI = min(x):.5:max(x);
YI = min(y):.5:max(y);
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