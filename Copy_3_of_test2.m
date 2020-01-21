% PLOTTING

clear; close all;
data1 = readtable('batch.csv');
data2 = readtable('20171124 MagCoord3axisData.csv');
% ---------------------------------------------
% idx = 1:195;
idx = 1:length(data1.y);
yMin = min(data1.y(idx));
yMax = max(data1.y(idx));
% lM = data1(data1.y>=yMin,:);
lM = table2array(data1(data1.y==21.1,:));
sortedlM = sortrows(lM,1);
% ---------------------------------------------
tM = [data2.mag_x,data2.mag_y,data2.mag_z];
% agl = -pi/2;
% agl = 0.1363;
% R = [cos(agl) -sin(agl) 0 ;
%     sin(agl) cos(agl) 0;
%     0 0 1];
% rotlM = (R*tM')';


% x = 5; % data's optimal value
x = -pi/2;
R = [cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1];
tM = (R*tM')'-25;


x1 = lM(:,1);
x2 = data2.coord_x;

mx1 = lM(:,3);
my1 = lM(:,4);
mz1 = lM(:,5);

mx2 = tM(:,1);
my2 = tM(:,2);
mz2 = tM(:,3);

subplot(411)
plot(x1,mx1,'s-', x2,mx2,'s-')
legend('Learning data', 'Test data')
title('Comparison of x magnet data')
subplot(412)
plot(x1,my1,'s-', x2,my2,'s-')
legend('Learning data', 'Test data')
title('Comparison of y magnet data')
subplot(413)
plot(x1,mz1,'s-', x2,mz2,'s-')
legend('Learning data', 'Test data')
title('Comparison of z magnet data')
subplot(414)
plot(x1,vecnorm([mx1 my1 mz1],2,2),'s-', x2,vecnorm([mx2 my2 mz2],2,2),'s-')
legend('Learning data', 'Test data')
title('Comparison of L2Norm magnet data')

set(gcf,'units','points','position',[400,200,1000,800])
sdf(gcf','sj')

%%
radian_x = linspace(-pi,pi,63);
%% --------------------------------------------- Calculate Data for 2D plot
figure
R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]),(radian_x)','UniformOutput',false);
rotatedMag = cell2mat(cellfun(@(x)((x*tM(1,:)')'),R,'UniformOutput',false));
mag_dist = pdist([lM(1,3:5);rotatedMag]);
mag_dist = mag_dist(1:length(R));

plot(radian_x,mag_dist)
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'}) 
set(gcf,'units','points','position',[800,500,1000,500])
sdf(gcf','sj')
% [~,I] = min(mag_dist);
% radian_x(I)
%% --------------------------------------------- Calculate Data for 3D plot
[X,Y] = meshgrid(x2, radian_x);
% [X,Y] = meshgrid(radian_x,x2);
Z = arrayfun(@(x,y) getDistance(x,[data2.coord_x,data2.coord_y,tM],lM,y),X,Y);
% R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]),Y,'UniformOutput',false);
% rotatedMag = cell2mat(cellfun(@(x,y)((y*tM(x,:)')'),mat2cell(X),R,'UniformOutput',false));
%% DRAWING 3D surf
figure
%%
norm_Z = smoothn(Z,'robust');
s = surf(X,Y,norm_Z,'FaceAlpha',0.5);
view(20.9000,57.2000)
% view(.5,90)
set(gca,'YTick',-pi:pi/2:pi) 
set(gca,'YTickLabel',{'-pi','-pi/2',0,'pi/2','pi'}) 
set(gcf,'units','points','position',[500,500,1200,800])
sdf(gcf','sj')
% s.EdgeColor = 'none';
% colormap(flipud(spring))
% colormap hot
% caxis([20 200])
colorbar
% colormap(parula(5))
%%

%%
% VISUAL EFFECT
% light;
% lighting phong;
% camlight('left');
% shading interp;
%%
function d = getDistance(x,tdata,ldata,rot_angle)
l_x = interp1(ldata(:,1),ldata(:,3),x);
l_y = interp1(ldata(:,1),ldata(:,4),x);
l_z = interp1(ldata(:,1),ldata(:,5),x);
t_x = interp1(tdata(:,1),tdata(:,3),x);
t_y = interp1(tdata(:,1),tdata(:,4),x);
t_z = interp1(tdata(:,1),tdata(:,5),x);
rotMat = @(x) ([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]); 
d = pdist([l_x,l_y,l_z;(rotMat(rot_angle)*[t_x;t_y;t_z])']);
end


