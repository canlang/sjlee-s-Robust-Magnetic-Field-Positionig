mag1 = readtable('db_20171020_15_26_24/magnetic.csv');
mag2 = readtable('db_20171020_15_27_55/magnetic.csv');
mag3 = readtable('db_20171020_15_29_21/magnetic.csv');
mag4 = readtable('db_20171020_15_30_46/magnetic.csv');

% acc = readtable('14_55_49/accelerometer.csv');

%%
subplot(811)
plot(mag1.time, mag1.x,mag1.time, mag1.y,mag1.time, mag1.z)
subplot(812)
plot(mag1.time, sqrt(sum([mag1.x mag1.y mag1.z].^2,2)))

subplot(813)
plot(mag2.time, mag2.x,mag2.time, mag2.y,mag2.time, mag2.z)
subplot(814)
plot(mag2.time, sqrt(sum([mag2.x mag2.y mag2.z].^2,2)))

subplot(815)
plot(mag3.time, mag3.x,mag3.time, mag3.y,mag3.time, mag3.z)
subplot(816)
plot(mag3.time, sqrt(sum([mag3.x mag3.y mag3.z].^2,2)))

subplot(817)
plot(mag4.time, mag4.x,mag4.time, mag4.y,mag4.time, mag4.z)
subplot(818)
plot(mag4.time, sqrt(sum([mag4.x mag4.y mag4.z].^2,2)))



% subplot(312)
% plot(acc.time, acc.x,acc.time, acc.y,acc.time, acc.z)