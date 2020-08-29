clear;close all;clc;

% ------------------------------------------------------------------------
% path_dir = 'rawdata';
% path_dir = 'rawdata/test-ki-huawei';
% path_dir = 'rawdata/test-KI-1F-MATE20pro';

% target_rawdata_paths = getNameFolds(path_dir);
% rawdata = load_rawdata(fullfile(path_dir,target_rawdata_paths{1})); 

% ------------------------------------------------------------------------
target_rawdata_paths = getNameFolds('rawdata');
% rawdata = load_rawdata(fullfile('rawdata',target_rawdata_paths{63})); 
% rawdata = load_rawdata(fullfile('rawdata',target_rawdata_paths{61})); 
% rawdata = load_rawdata(fullfile('rawdata',target_rawdata_paths{43})); 
% rawdata = load_rawdata(fullfile('rawdata',target_rawdata_paths{16})); % for ieee accesS?
% rawdata = load_rawdata(fullfile('rawdata',target_rawdata_paths{86}));
rawdata = load_rawdata(fullfile('rawdata',target_rawdata_paths{91}));
% index 9 also available
% index 36 good sample (tracking visual)

% rawdata = load_rawdata('181217_170914_656_N1_긴복도_동쪽방향');
% rawdata = load_rawdata('181217_171014_406_N1_긴복도_서쪽방향');
% rawdata = load_rawdata('rawdata/190109_133853_526_N1_연구실시작종료');
% rawdata = load_rawdata('rawdata/190110_101202_128_N1_연구실시작종료');
% rawdata = load_rawdata('rawdata/190110_102026_029_N1_연구실시작종료');
% rawdata = load_rawdata('rawdata/200714_132017_914_doan_straight_return');
% rawdata = load_rawdata('rawdata/200715_160602_093_n1_r_turn_return');
% rawdata = load_rawdata('rawdata/200717_143220_867_n1_fullcover');

% rawdata = load_rawdata('input_rawdata/190409_172013_390_N1_upper_corridor_rover');

%% resample
% (1)
% raw_acc = rawdata.acc;
% raw_gyr = rawdata.gyr;
% acc_time = rawdata.acc_time;
% gyr_time = rawdata.gyr_time;
% acc_mag = rawdata.acc_norm;
% T_acc = timetable(seconds(raw_acc(:,2)/1e9),raw_acc(:,3:5),acc_mag);
% T_gyr = timetable(seconds(raw_gyr(:,2)/1e9),raw_gyr(:,3:5));
% T_acc = sortrows(T_acc);
% T_gyr = sortrows(T_gyr);
% 
% TT = synchronize(T_acc,T_gyr,'regular','linear','TimeStep',seconds(2e-2));
% % TT = synchronize(T_acc,T_gyr,'commonrange','linear','TimeStep',seconds(2e-2));
% % TT = synchronize(T_acc,T_gyr,'intersection','TimeStep',seconds(2e-2));
% % TT = synchronize(T_acc,T_gyr,'commonrange','linear');
% 
% % TT = synchronize(T_acc,T_gyr);
% % TT = synchronize(T_acc,T_gyr,'regular','nearest','TimeStep',seconds(2e-2));
% 
% % TT = synchronize(T_acc,T_gyr,'commonrange','SampleRate',50,'method','nearest');
% % TT = synchronize(T_acc,T_gyr,'intersection','mean');
% 
% Accelerometer = TT.Var1_T_acc;
% Gyroscope = TT.Var1_T_gyr*180/pi;
% acc_mag = TT.acc_mag;
% 
% time = seconds(TT.Time(:)-(TT.Time(1)));
% rate = median(diff(time)); % cal sample rate
% 
% % acc = resample(acc_mag,time,20);
% % [Accelerometer,time] = resample(acc_mag, raw_acc(:,2)/1e9,20);

% (2) 2020.Aug.29
sample_rate = 2e-2;                                     % reampling rate
res_data = resample_rawdata2(rawdata, sample_rate);
acc_mag = res_data.acc_norm;
rate = sample_rate;
time = seconds(res_data.Time(:)-(res_data.Time(1)));
Accelerometer = res_data.Accelerometer;
Gyroscope = res_data.Gyroscope*180/pi;

%% find step point (step) and time labeling
% threshold should be tuned experimentally to match a person's level
minPeakHeight = std(acc_mag);

[pks,locs] = findpeaks(acc_mag,'MinPeakDistance',...
    .3/rate,'MinPeakHeight',minPeakHeight);   % .3s 이내의 피크는 무시

%%  
addpath(genpath('madgwick_algorithm_matlab'));
AHRS = MadgwickAHRS('SamplePeriod', rate, 'Beta', 0.1); % sample rate: 2e-2
% AHRS = MahonyAHRS('SamplePeriod', rate, 'Kp', 0.1); % sample rate: 2e-2

quaternion = zeros(length(time), 4);
for t = 1:length(time)
    AHRS.Update(Gyroscope(t,:) * (pi/180), Accelerometer(t,:), Accelerometer(t,:));	% gyroscope units must be radians
    quaternion(t, :) = AHRS.Quaternion;
end

euler = quatern2euler(quaternConj(quaternion)) * (180/pi);	% use conjugate for sensor frame relative to Earth and convert to degrees.

%%
estloc = zeros(length(locs),2);
% phi,theta,psi: roll pitch yaw (aerospace sequence?)
for i=1:length(locs)
% for i=1:length(estloc)
    t_i = locs(i);
    yaw = euler(t_i,3);
%     yaw = euler(t_i,3)+90;
    step = 0.7;
    trM = [step*(1*cosd(yaw) - 0*sind(yaw));step*(1*sind(yaw) + 0*cosd(yaw))];
%     trM = [step*(0*cosd(yaw) - 1*sind(yaw));step*(0*sind(yaw) + 1*cosd(yaw))];
    if all(estloc == 0)
        estloc(i,:) = (trM+[0;0])';
    else
        estloc(i,:) = (trM+estloc(i-1,:)')';
    end
end

%% draw v.02
close all
figure
subplot(211)
yaw = unwrap((euler(:,3))*pi/180);
% s_yaw = smoothdata(yaw,'movmedian',0);
% plot(time, deg2rad(yaw))
plot(time, (yaw))
hold on
% yline(2.5);
yl = yline(0,'--','y = 0','LineWidth',3);
% yl = yline(2*pi,'g--','LineWidth',3);
set(gca,'YTick',0:pi/2:2*pi) 
set(gca,'YTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'}) 
ylim([-pi/2, 2*pi])
xlabel('Time (sec)');ylabel('\psi (rad)')
% title('time to yaw (rad)')

% subplot(212)
% plot(time, stdfilt(yaw),'.')
% xlabel('time (sec)');ylabel('\sigma of \psi (rad)')
% % title('time to std (yaw)')
% % plot((euler(:,3)))
% set(gcf,'units','points','position',[500,500,800,500])
% sdf(gcf,'sj2')

subplot(212)
% figure
% subplot(3,2,2:2:4)
plot(estloc(:,1),estloc(:,2),'xm-','MarkerSize',8)
xlabel('x (m)');
ylabel('y (m)');
% ylim([-10 50])
axis image
% axis 'auto'
grid on
% grid minor
ax = gca;
% set(gcf,'units','points','position',[800, 100, 641, 288])
set(gcf,'units','points','position',[1096         499         659         673])
% set(gcf,'units','points','position',[1300,500,800,600])
sdf(gcf,'sj2')

print -clipboard -dbitmap
% print -depsc2 eps/18times_repeat_circle_pdr.eps
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