clear; close all;
data1 = readtable('batch.csv');
data2 = readtable('20171124 MagCoord3axisData.csv');
%% ---------------------------------------------
% idx = 1:195;
idx = 1:length(data1.y);
yMin = min(data1.y(idx));
yMax = max(data1.y(idx));
% lM = data1(data1.y>=yMin,:);
lM = table2array(data1(data1.y==21.1,:));
sortedlM = sortrows(lM,1);
%%
% plot(data1.x,data1.y,'x')
% hold on
% plot(lM(:,1),lM(:,2),'x')
% hold off
% axis equal
subplot(311)
x = 1:length(lM);
% plot(x,lM(:,3))
plot(x,lM(:,3),x,lM(:,4),x,lM(:,5))
title('Learning data : magnet x')
% hold on
% plot(data2.coord_x,data2.coord_y,'o')
% hold off
% axis equal
% subplot(311)
% plot(lM.magnet_x)
%% ---------------------------------------------
tM = [data2.mag_x,data2.mag_y,data2.mag_z];
agl = -pi/2;
R = [cos(agl) -sin(agl) 0 ;
    sin(agl) cos(agl) 0;
    0 0 1];
rotlM = (R*tM')';

x = -pi/2;
R = [cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1];
tM = (R*tM')';

subplot(312)
x = 1:63;
plot(x, data2.mag_x,x, data2.mag_y,x, data2.mag_z)
title('Raw test data magnet')
subplot(313)
plot(x,tM(:,1),x,tM(:,2),x,tM(:,3))
title('Rotated test data magnet')
%% ---------------------------------------------
