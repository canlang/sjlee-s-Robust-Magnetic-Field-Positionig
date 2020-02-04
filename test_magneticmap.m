clear; close all;
% addpath('xyz file operations')

%%
% data1 = readtable('batch.csv');
% lM = [data1.magnet_x,data1.magnet_y,data1.magnet_z];
% x = data1.x;
% y = data1.y;
% z = vecnorm(lM,2,2);
%%
% data1 = readtable('sw_maps/dataset_n1_2f_result.csv');
% lM = [data1.mag_x,data1.mag_y,data1.mag_z];
% loc_x = data1.loc_x;
% loc_y = data1.loc_y;
% z = vecnorm(lM,2,2);
%%
% data1 = readtable('sw_maps/dataset_ki_1f.csv');
% data1 = rmmissing(data1);
% lM = [data1.mag_x,data1.mag_y,data1.mag_z];
% x = data1.loc_x;
% y = data1.loc_y;
% z = vecnorm(lM,2,2);
%%
load('mats/magmap-n1-2f-0.6p.mat');
x = map(:,1);
y = map(:,2);
lM = map(:,3:5);
% lM = [data1.mag_x,data1.mag_y,data1.mag_z];
% x = data1.loc_x;
% y = data1.loc_y;
z = vecnorm(lM,2,2);
%%
[X,Y,Z] = xyz2grid(x,y,z);
% plot(loc_x,loc_y,'.')
%%
% XYZ = [x,y,z];
% % Get coordinate vectors
% ux = unique(XYZ(:,1)) ;
% uy = unique(XYZ(:,2)) ;
% % dimensions of the data
% nx = length(ux) ; 
% ny = length(uy) ;
% % Frame matrix of grid 
% D = reshape(XYZ(:,3),[ny,nx]) ;
% % flip  matrix to adjust for plot
% H = flipud(H) ;
% % Transpose the matrix 
% H = H' ;  % Check if is required
% axis image
% axis xy
%% imagesc version

ax1 = subplot(211);
h = imagesc(X(1,:),Y(:,1),Z);

xlabel('x (m)')
ylabel('y (m)')
% axis image
% axis equal
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

XI = unique(x);
YI = unique(y);
[X,Y] = meshgrid(XI,YI);
Z = griddata(x,y,vecnorm(lM,2,2),X,Y);
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