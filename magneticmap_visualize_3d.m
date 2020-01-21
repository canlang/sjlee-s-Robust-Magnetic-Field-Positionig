clear; close all;

% addpath('xyz file operations')

% data1 = readtable('batch.csv');
% lM = [data1.magnet_x,data1.magnet_y,data1.magnet_z];

data1 = readtable('sw_maps/dataset_ki_1f.csv');
data1 = rmmissing(data1);
lM = [data1.mag_x,data1.mag_y,data1.mag_z];


% DATA SPLIT
x = data1.loc_x;
y = data1.loc_y;
% z{1} = vecnorm(lM,2,2);    % L2 norm
z{1} = data1.mag_x;     % x component
z{2} = data1.mag_y;     % y component
z{3} = data1.mag_z;     % z component
hold on
for k=1:3
    [X,Y,Z] = interpolation_by_alphashape(x,y,z{k},.3);
%     Z1 = Z;
%     Z1(~isnan(Z)) = 1;  
    surface('XData',X, 'YData',Y, 'ZData',Z.*.03+k*15, ...
            'CData',Z, 'CDataMapping','scaled', ...
            'EdgeColor','none', 'FaceColor','texturemap','FaceAlpha',1)
    colormap hsv
    freezeColors
%     annotation('textarrow',[30,40],[k*15+1,k*15],'String','y = x ')
%     annotation(figure1,'textarrow',[0.848576779026218 0.801910112359552],...
%     [0.655341463414635 0.629674796747968],
end

view(3), axis tight image
view(25,34)
set(gca,'xtick',[])
set(gca,'ytick',[])
set(gca,'ztick',[])


I = imread('N1-7F.png','BackgroundColor',[1 1 1]);
I = flipud(I);
[X,Y] = meshgrid(1:size(I,2), 1:size(I,1));
Z = ones(size(I,1),size(I,2));


xWorldLimits = [-1 1650/20];
yWorldLimits = [-1 660/20];
RA = imref2d(size(I),xWorldLimits,yWorldLimits);
surface('XData',X/20-0.1, 'YData',Y/20-0.5, 'ZData',Z.*0, ...
        'CData',I, 'CDataMapping','direct', ...
        'EdgeColor','none', 'FaceColor','texturemap')

set(gcf,'units','points','position',[500,500,1200,600])
sdf(gcf,'sj2')

% hold on

% surface('XData',X-0.5, 'YData',Y-0.5, 'ZData',Z.*1, ...
%         'CData',I, 'CDataMapping','direct', ...
%         'EdgeColor','none', 'FaceColor','texturemap')
    
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
return
%%
h = contourf(X,Y,Z);
hcb = colorbar;
ylabel(hcb, 'L2 norm of Magnetic Field')
% ylabel(hcb, 'X compoment of Magnetic Field')
% ylabel(hcb, 'Y compoment of Magnetic Field')
% ylabel(hcb, 'Z compoment of Magnetic Field')

xlabel('x (m)')
ylabel('y (m)')
axis image
axis xy
%% ----------------------------------------------------
set(gcf,'units','points','position',[500,500,1200,600])
sdf(gcf,'sj2')
%% ----------------------------------------------------
% saveas(gcf,'eps/interp_norm.eps','epsc');
% saveas(gcf,'eps/interp_x.eps','epsc');
% saveas(gcf,'eps/interp_y.eps','epsc');
% saveas(gcf,'eps/interp_z.eps','epsc');

%% INTERPOLATION by alphashape
function [X,Y,Z] = interpolation_by_alphashape(x,y,z,gt)
XI = min(x):gt:max(x);
YI = min(y):gt:max(y);
[X,Y] = meshgrid(XI,YI);

% APPLY ALPHA SHAPE 
shp = alphaShape(x,y);
% shp.Alpha = 2;            % Alpha parameter: 
in = inShape(shp,X,Y);
xg = X(in);
yg = Y(in);
% DATA INPUT ONLY INSIDE OF ALPHA SHAPE
zg = griddata(x,y,z,xg,yg);         % 2. griddata() : INTERPOLATION
[X,Y,Z] = xyz2grid(xg,yg,zg);
end