clear;

mag1 = readtable('db_20171020_15_47_05/magnetic.csv');
mag2 = readtable('db_20171020_15_48_22/magnetic.csv');
mag3 = readtable('db_20171020_15_49_38/magnetic.csv');
mag4 = readtable('db_20171020_15_51_04/magnetic.csv');

% acc = readtable('14_55_49/accelerometer.csv');

rotZrad = deg2rad(180);
rotZ = [cos(rotZrad), -sin(rotZrad), 0;
    sin(rotZrad), cos(rotZrad), 0;
    0, 0, 1];


%%
subplot(411)
plot(mag1.time, mag1.x,mag1.time, mag1.y,mag1.time, mag1.z)
subplot(412)
plot(mag1.time, sqrt(sum([mag1.x mag1.y mag1.z].^2,2)))

subplot(413)
plot(mag3.time, -mag3.x,mag3.time, -mag3.y,mag3.time, mag3.z)
subplot(414)
plot(mag3.time, sqrt(sum([mag3.x mag3.y mag3.z].^2,2)))



% subplot(312)
% plot(acc.time, acc.x,acc.time, acc.y,acc.time, acc.z)

%% <- right to left (way)
subplot(311)
plot(1:length(mag1.time), mag1.x,1:length(mag3.time), -mag3.x)
subplot(312)
plot(1:length(mag1.time), mag1.y,1:length(mag3.time), -mag3.y)
subplot(313)
plot(1:length(mag1.time), mag1.z,1:length(mag3.time), mag3.z)

%% -> left to right (way)
% rotedMag = rotZ*[mag4.x'; mag4.y'; mag4.z'];
% rotedMag = rotedMag';
% subplot(311)
% plot(1:length(mag2.time), mag2.x,1:length(mag4.time), rotedMag(:,1))
% subplot(312)
% plot(1:length(mag2.time), mag2.y,1:length(mag4.time), rotedMag(:,2))
% subplot(313)
% plot(1:length(mag2.time), mag2.z,1:length(mag4.time), rotedMag(:,3))

% subplot(311)
% plot(1:length(mag2.time), mag2.x,1:length(mag4.time), -mag4.x)
% subplot(312)
% plot(1:length(mag2.time), mag2.y,1:length(mag4.time), -mag4.y)
% subplot(313)
% plot(1:length(mag2.time), mag2.z,1:length(mag4.time), mag4.z)