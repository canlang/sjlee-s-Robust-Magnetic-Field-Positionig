clearvars; close all;
% rawdata = load_rawdata('rawdata/200714_132017_914_doan_straight_return');
% rawdata = load_rawdata('rawdata/190109_133853_526_N1_연구실시작종료');
% rawdata = load_rawdata('rawdata/200715_160602_093_n1_r_turn_return');
rawdata = load_rawdata('rawdata/200717_143220_867_n1_fullcover');

% rawdata = load_rawdata('input_rawdata/190409_172013_390_N1_upper_corridor_rover');

% --------------------------------------------------------------
sample_rate = 2e-2;                                     % reampling rate
res_data = resample_rawdata2(rawdata, sample_rate);     % resample data

% --------------------------------------------------------------
min_peak_thr = std(res_data.acc_norm);
[~,step_idx] = findpeaks(res_data.acc_norm,'MinPeakDistance',.3/sample_rate,...
    'MinPeakHeight',min_peak_thr);
% visualize for check the step state
% findpeaks(res_data.acc_norm,'MinPeakDistance',.3/res_rate,...   %
%     'MinPeakHeight',min_peak_thr)
% --------------------------------------------------------------
addpath(genpath('madgwick_algorithm_matlab'));                  % generate rotation matrix to (transform?)
AHRS = MadgwickAHRS('SamplePeriod', sample_rate, 'Beta', 0.1);  % sample rate: 2e-2
quaternion = zeros(length(res_data.Time), 4);
for t = 1:length(res_data.Time)
    AHRS.Update(res_data.Gyroscope(t,:), ...
        res_data.Accelerometer(t,:), ...
        res_data.Accelerometer(t,:));	% gyroscope units must be radians
    quaternion(t, :) = AHRS.Quaternion;
end
% --------------------------------------------------------------

% R_matrix = quatern2rotMat(quaternion(step_idx,:));
R_matrix = quatern2rotMat(quaternConj(quaternion(step_idx,:)));
euler = quatern2euler(quaternConj(quaternion)) * (180/pi);
% euler = quatern2euler(quaternion) * (180/pi);

input_mag = res_data.Magnetometer(step_idx,:);

output_mag = zeros(size(input_mag));
for i=1:length(input_mag)
    output_mag(i,:) = (R_matrix(:,:,i)*input_mag(i,:)')';
end

% 1. not optimized
% estloc = getTrajectory(euler,step_idx);
% 2. optimized
load('trajectory/gt_traj_full')
[estloc, c_lambda, ~] = getOptTrajectory(euler,step_idx,gt);
% cum_lambda = zeros(size(lambda));
for i=1:length(c_lambda)
%     cum_lambda(i:end,1) = cum_lambda(i:end,1) + lambda(i);
    el = c_lambda(i);
    output_mag(i,:) = ([cosd(el) -sind(el) 0;...
                        sind(el)  cosd(el) 0;...
                        0         0        1]*output_mag(i,:)')';
end

map = [estloc output_mag];

mat_filename = fullfile('mats',['magmap-n1-7f-','step','w.mat']);
if ~exist(mat_filename,'file')
    save(mat_filename,'map')
end

%% plotting and visualize trajectory
subplot(2,1,1)
sc = 2;     % selected component
time_x = datetime(res_data.Time,'ConvertFrom','posixtime','TimeZone','Asia/Seoul');
% plot(raw_data.mag(:,3:5))
plot(time_x, res_data.Magnetometer(:,sc),'-')
hold on
plot(time_x(step_idx), input_mag(:,sc), '.');
hold off
legend('rawdata mag_x','step')

subplot(2,1,2)
plot(time_x, res_data.Magnetometer(:,sc),'-')
hold on
plot(time_x(step_idx), output_mag(:,sc), '.');
hold off
legend('rawdata mag_x','R(mag_x)')

set(gcf,'units','points','position',[300,500,1200,600])
% sdf(gcf,'sj2')

figure
% subplot(3,2,2:2:4)
plot(estloc(:,1),estloc(:,2),'xm-','MarkerSize',8)
xlabel('x (m)');
ylabel('y (m)');
% ylim([-10 50])
axis image
% axis 'auto'
grid on
grid minor
