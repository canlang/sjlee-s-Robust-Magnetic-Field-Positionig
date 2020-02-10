clearvars;close all;

tidx = 1;
files = dir('rawdata/test-n1-7f-iphone/*.csv');
file_folder = files(tidx).folder;
file_name = files(tidx).name;
path = fullfile(file_folder,file_name);

% TT = readtimetable(path);
% M = readmatrix(path);
T = readtable(path,'Delimiter',',','PreserveVariableNames',true);
% T(1,1)
dt = datetime(T{:,1},'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS +0900');
pt = posixtime(dt);

% T{:,40:42}+T{:,48:50}

acc = [T.("accelerometerAccelerationX(G)"),T.("accelerometerAccelerationY(G)"),T.("accelerometerAccelerationZ(G)")]...
        +[T.("motionGravityX(G)"),T.("motionGravityY(G)"),T.("motionGravityZ(G)")];

subplot(411)
mag = [T.("motionMagneticFieldX(µT)"),T.("motionMagneticFieldY(µT)"),T.("motionMagneticFieldZ(µT)")];
plot(mag)
subplot(412)
mag = [T.("locationHeadingX(µT)"),T.("locationHeadingY(µT)"),T.("locationHeadingZ(µT)")];
plot(mag)
subplot(413)
mag1 = readtable('db_20171020_15_47_05/magnetic.csv');
plot(mag1{:,3:5})
subplot(414)
mag3 = readtable('db_20171020_15_49_38/magnetic.csv');
plot(mag3{:,3:5})
