clear; close all;
data1 = readtable('batch.csv');
data2 = readtable('20171124 MagCoord3axisData.csv');
data3 = readtable('shatie_testdata_N1_20180712/00.csv');
data4 = readtable('sjlee_testdata_N1_20180723/sjlee_testdata_N1_20180723.csv');


tM = [data2.mag_x,data2.mag_y,data2.mag_z];
% rotate test data 
x = -pi/2;
R = [cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1];
tM = (R*tM')'-25;

subplot(311)
plot(data2.coord_x,tM(:,1),data4.location_x, data4.magnet_x)
subplot(312)
plot(data3.magnet_x)
subplot(313)
plot(data4.location_x, data4.magnet_x)