clear; close all;
addpath('xyz file operations')

data1 = readtable('batch.csv');

lM = [data1.magnet_x,data1.magnet_y,data1.magnet_z];
x = data1.x;
y = data1.y;
z = vecnorm(lM,2,2);
[X,Y,Z] = xyz2grid(x,y,z);
%% imagesc version
ax1 = subplot(211);
h = imagesc(X(1,:),Y(:,1),Z);

xlabel('x (m)')
ylabel('y (m)')
axis image
axis xy
set(h,'alphadata',~isnan(Z))
colormap(ax1,cool)
colorbar
%% heatmap version
% h = heatmap(Z,'GridVisible','off','MissingDataColor',[1 1 1]);
% heatmap(Z,[],[],[],'ColorMap', @cool, 'NaNColor', [0 0 0], 'colorbar', true);
% set(h,'alphadata',~isnan(Z))

%% interpolation
ax2 = subplot(212);
XI = unique(data1.x);
YI = unique(data1.y);
[X,Y] = meshgrid(XI,YI);
Z = griddata(data1.x,data1.y,vecnorm(lM,2,2),X,Y);
h = imagesc(XI,YI,Z);
xlabel('x (m)')
ylabel('y (m)')
axis image
axis xy
set(h,'alphadata',~isnan(Z))
colormap(ax2,cool)
colorbar
% F = scatteredInterpolant(data1.magnet_x,data1.magnet_y,data1.magnet_z,'nearest');
% vq = F(X,Y);

%%
set(gcf,'units','points','position',[100,500,1200,600])
sdf(gcf,'sj2')

%%
% figure
% surf(X,Y,Z)

% heatmap(XI,YI,vq,'GridVisible','off')
% axis equal
%%
% xvector = x;
% yvector = y;
% numX = numel(xvector);  
% numY = numel(yvector);
% yvector = repmat(yvector(:),numX,1);
% xvector = repmat(xvector   ,numY,1);
% XY = [xvector(:) yvector];
% Z=(1:1:length(XY))';

% XYZ=[x,y,z];
% M=accumarray(XYZ(:,[2 1]), XYZ(:,3));
% [X,Y] = ind2sub(size(M), 1:numel(M));
% XYZ = [X(:), Y(:), M(:)];
% %%
% heatmap(XYZ,'GridVisible','off')