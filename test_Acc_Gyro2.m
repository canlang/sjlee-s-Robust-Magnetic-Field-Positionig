clear;close all;clc;
addpath('xyz file operations')
addpath('Multiprod_2009')


target_rawdata_paths = getNameFolds('input_rawdata');
rawdata = load_rawdata(fullfile('input_rawdata',target_rawdata_paths{23}));


%% resample for synchronize
rate = 2e-2;
processed_data = resample_rawdata(rawdata,rate);


%% find step point (step) and time labeling
% threshold should be tuned experimentally to match a person's level
minPeakHeight = std(processed_data.acc_norm);

[pks,locs] = findpeaks(processed_data.acc_norm,'MinPeakDistance',...
    .3/rate,'MinPeakHeight',minPeakHeight);   % .3s 이내의 피크는 무시

%% rebuild trajectory
addpath(genpath('madgwick_algorithm_matlab'));
AHRS = MadgwickAHRS('SamplePeriod', rate, 'Beta', 0.1); % sample rate: 2e-2

quaternion = zeros(length(processed_data.Time), 4);
for t = 1:length(processed_data.Time)
    AHRS.Update(processed_data.Gyroscope(t,:), processed_data.Accelerometer(t,:),...
        processed_data.Accelerometer(t,:));	% gyroscope units must be radians
    quaternion(t, :) = AHRS.Quaternion;
end

euler = quatern2euler(quaternConj(quaternion)) * (180/pi);	% use conjugate for sensor frame relative to Earth and convert to degrees.

%% labeling location 
estloc = zeros(length(locs),2);
% phi,theta,psi: roll pitch yaw (aerospace notation sequence?)
for i=1:length(locs)
% for i=1:length(estloc)
    t_i = locs(i);
    yaw = euler(t_i,3)+0;
    step = 0.7;
    trM = [step*(1*cosd(yaw) - 0*sind(yaw));step*(1*sind(yaw) + 0*cosd(yaw))];
%     trM = [step*(0*cosd(yaw) - 1*sind(yaw));step*(0*sind(yaw) + 1*cosd(yaw))];
    if all(estloc == 0)
        estloc(i,:) = (trM+[0;0])';
    else
        estloc(i,:) = (trM+estloc(i-1,:)')';
    end
end

%% temporal
% theta = unwrap(deg2rad(euler(locs,3)));
theta = deg2rad(euler(locs,3))+0*(pi/180);
r = 1*ones(size(theta));
% vecnorm(euler(locs,:),2,2)

[u,v] = pol2cart(theta,r);
% feather(u,v)
quiver(estloc(:,1),estloc(:,2),u,v,.4)
axis image
grid on

% set(gcf,'units','points','position',[500,500,800,1000])
set(gcf,'units','points','position',[500,500,1200,800])
sdf(gcf,'sj2')

print -clipboard -dbitmap
return
%% draw v.02
close all
figure
subplot(211)
% yaw = unwrap((euler(:,3)));  % effective???
yaw = (euler(:,3));
% s_yaw = smoothdata(yaw,'movmedian',0);
plot(processed_data.Time, deg2rad(yaw))
xlabel('time (sec)');ylabel('\psi (rad)')
% title('time to yaw (rad)')
subplot(212)
plot(processed_data.Time, stdfilt(deg2rad(yaw)),'.')
xlabel('time (sec)');ylabel('\sigma of \psi (rad)')
% title('time to std (yaw)')
% plot((euler(:,3)))
set(gcf,'units','points','position',[500,500,800,500])
sdf(gcf,'sj2')

figure
% subplot(3,2,2:2:4)
plot(estloc(:,1),estloc(:,2),'xr-','MarkerSize',8)
xlabel('m');
ylabel('m');
% ylim([-10 50])
axis image
% axis 'auto'
grid on

set(gcf,'units','points','position',[500,500,800,1000])
sdf(gcf,'sj2')

print -clipboard -dbitmap
return
%% draw v.01
figure
subplot(325)
plot(acc_time, raw_acc(:,3:5))

% subplot(312)
% plot(gyr_time, raw_gyr(:,3:5))

subplot(321)
yaw = unwrap((euler(:,3)));
% s_yaw = smoothdata(yaw,'movmedian',0);
plot(time, deg2rad(yaw))
% xlabel('time');ylabel('yaw')
title('time to yaw (rad)')
subplot(323)
plot(time, stdfilt(deg2rad(yaw)))
title('time to std (yaw)')
% plot((euler(:,3)))

subplot(3, 2, 2:2:6)
plot(estloc(:,1),estloc(:,2),'xr-','MarkerSize',8)
% ylim([-10 50])
axis image
% axis 'auto'
grid on

set(gcf,'units','points','position',[500,500,1200,800])
sdf(gcf,'sj2')

print -clipboard -dbitmap
%% local functions
function nameFolds = getNameFolds(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
% nameFolds(~cellfun(@isempty,regexp(nameFolds, '\d{6}'))) = [];
end