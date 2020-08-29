clear; close all;
% batch script to create eps files

addpath('xyz file operations')


site_name = 'KI-1F';
% site_name = "N1-7F";
A = imread(sprintf('map/%s.png',site_name),'BackgroundColor',[1 1 1]);

switch site_name
    
    case 'KI-1F'
        map = magmap_construction('mats',site_name,.1);
        x = map(:,1);
        y = map(:,2);
        lM = map(:,3:5);
        
        scale = 0.05;
        xWorldLimits = [0 size(A,2)*scale];
        yWorldLimits = [0 size(A,1)*scale];
    case 'N1-7F'
        data1 = readtable('batch.csv');
        x = data1.x;
        y = data1.y;
        lM = [data1.magnet_x,data1.magnet_y,data1.magnet_z];
        
        xWorldLimits = [-1 1650/20];
        yWorldLimits = [-1 660/20];
end
data_seq = ["x","y","z","norm"];
for i=1:4
    % DATA SPLIT
    if i==4
        z = vecnorm(lM,2,2);    % L2 norm
    else
        z = lM(:,i);     % x component
        % z = lM(:,2);     % y component
        % z = lM(:,3);     % z component
    end
    

    figure
    RA = imref2d(size(A),xWorldLimits,yWorldLimits);
    imshow(flipud(A),RA);
%     axis image
% N1
%     xlim([5   83])
%     ylim([5.000-3   25.000+3])    
% KI
    xlim([3.0000-15   50.5000+15])
    ylim([7.5000-3   45.5000+3])
    axis xy
%     grid on

    hold on

    %%
    if site_name == "KI-1F"
        [X,Y,Z] = interpolation_by_alphashape(x,y,z,0.5);
    else
        [X,Y,Z] = interpolation_by_alphashape(x,y,z);
    end


    % When alpha shape draw (need modify function in interpolation_by_alphashape)
    % set(gca,'XTick',[]);
    % set(gca,'YTick',[]);
    % sdf(gcf,'sj4')
    % print -depsc2 eps/alphashape.eps

    %%
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

    caxis([-100, 100]);
        
    hcb = colorbar;
    ylabel(hcb, [data_seq(i), ' of Magnetic Field (\muT)'])

    % ylabel(hcb, 'X compoment of Magnetic Field')
    % ylabel(hcb, 'Y compoment of Magnetic Field')
    % ylabel(hcb, 'Z compoment of Magnetic Field')


%% ----------------------------------------------------
    xlabel('x (m)')
    ylabel('y (m)')


% set(gcf,'units','points','position',[700,500,1000,350])


% set(gcf,'units','points','position',[500,500,1200,600])


%% ----------------------------------------------------
% saveas(gcf,'eps/interp_norm.eps','epsc');
% saveas(gcf,'eps/interp_x.eps','epsc');
% saveas(gcf,'eps/interp_y.eps','epsc');
% saveas(gcf,'eps/interp_z.eps','epsc');
% saveas(gcf,'eps/alphashape.eps','epsc');
% export_fig eps/alphashape.eps -depsc -m2
    sdf(gcf,'sj4')
%     tightfig
    saveas(gcf,sprintf('eps/%s-interp_%s.eps',site_name,data_seq(i)),'epsc');
end
%% INTERPOLATION by alphashape
function [X,Y,Z] = interpolation_by_alphashape(x,y,z,alpha)
XI = min(x):.5:max(x);
YI = min(y):.5:max(y);
[X,Y] = meshgrid(XI,YI);

% APPLY ALPHA SHAPE 
shp = alphaShape(x,y);
if nargin > 3   
    % shp.Alpha = .5;     % for KI building
    shp.Alpha = alpha;     % for KI building
end

% plot(shp,'EdgeColor','k')
% % plot(shp,'FaceAlpha',0.7,'EdgeAlpha',.7,'LineWidth',0.3)
% % plot(shp,'FaceColor','red','FaceAlpha',0.3,'LineWidth',0.3)
% % polyin = polyshape(shp.Points(:,1),shp.Points(:,2));
% % plot(polyin)
% 
% legend('alphashape')

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