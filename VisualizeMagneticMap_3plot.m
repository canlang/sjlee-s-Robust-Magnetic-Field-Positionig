clear; close all;
data1 = readtable('batch.csv');

lM = [data1.magnet_x,data1.magnet_y,data1.magnet_z];
for i=1:3
    % DATA SPLIT
    x = data1.x;
    y = data1.y;
    z = lM(:,i);     % x component

    %% INTERPOLATION
    XI = min(x):.5:max(x);
    YI = min(y):.5:max(y);
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

    A = imread('N1-7F.png','BackgroundColor',[1 1 1]);

    xWorldLimits = [-1 1650/20];
    yWorldLimits = [-1 660/20];
    RA = imref2d(size(A),xWorldLimits,yWorldLimits);
    subplot(3,1,i)
    imshow(flipud(A),RA);
    hold on

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
    h = contourf(X,Y,Z);
    hcb = colorbar;
    ylabel(hcb, 'L2 norm')
    xlabel('x (m)')
    ylabel('y (m)')
    axis image
    axis xy
end
%%
set(gcf,'units','points','position',[100,500,1200,1200])
sdf(gcf,'sj2')
% ----------------------------------------------------