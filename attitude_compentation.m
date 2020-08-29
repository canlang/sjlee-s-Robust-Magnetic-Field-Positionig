clearvars;close all;clc;
target_rawdata_paths = getNameFolds('rawdata');
% data_idxes = 87:89;
data_idxes = 92:-1:90;
% data_idxes = 95:-1:93;
atti = [0,45,85];

for i=1:3
    dataset_idx = data_idxes(i);
    dataname = target_rawdata_paths{dataset_idx};
    rawdata = load_rawdata(fullfile('rawdata',dataname));
    fprintf('dataset: %s\n',dataname)

    % if isequal(device_name, 'iphone')
    %     test_rawdata_paths = getNameFiles(test_path);
    %     rawdata = load_rawdata(fullfile(test_path,test_rawdata_paths{tr_idx}),'iPhone');    % iOS
    % else        
    %     test_rawdata_paths = getNameFolds(test_path);
    %     rawdata = load_rawdata(fullfile(test_path,test_rawdata_paths{tr_idx}));             % Android
    % end

    rate = 2e-2;
    processed_data = resample_rawdata2(rawdata,rate);  
    
    sidx = 1;
    processed_data.Accelerometer = processed_data.Accelerometer(sidx:end,:);
    processed_data.Gyroscope = processed_data.Gyroscope(sidx:end,:);
    processed_data.Magnetometer = processed_data.Magnetometer(sidx:end,:);
    processed_data.Time = processed_data.Time(sidx:end,:);
    
    addpath(genpath('madgwick_algorithm_matlab'));
    AHRS = MadgwickAHRS('SamplePeriod', rate, 'Beta', 0.1); % sample rate: 2e-2
%     AHRS = MahonyAHRS('SamplePeriod', rate, 'Kp', 0.1); % sample rate: 2e-2
    quaternion = zeros(length(processed_data.Time), 4);
    rot_mag = zeros(length(processed_data.Time), 3);
    for t = 1:length(processed_data.Time)
    %     AHRS.Update(processed_data.Gyroscope(t,:) * (pi/180), ... % This is
    %     bcause, in 'resample_rawdata' function make rad/s unit with Gyro.
        AHRS.Update(processed_data.Gyroscope(t,:), ...    
            processed_data.Accelerometer(t,:), ...
            processed_data.Accelerometer(t,:));	% gyroscope units must be radians
        quaternion(t, :) = AHRS.Quaternion;
%         rot_mag(t,:) = (quatern2rotMat(quaternConj(AHRS.Quaternion))*processed_data.Magnetometer(t,:)')';
        rot_mag(t,:) = quatrotate(quaternConj(AHRS.Quaternion),processed_data.Magnetometer(t,:));
    end
    
    subplot(3,4,1+4*(i-1))
    A = imread(sprintf('jpg/attitude_%d.jpg',atti(i)));    
    R = imrotate(A,-90);
    I = imcrop(R,[0 600 3024 2800]);
    imshow(I)
    title(sprintf('Attitude (θ) ≈ %d° ',atti(i)))
    
%     subplot(3,3,1+1*(i-1))  
    subplot(3,4,2+4*(i-1))
    mag_norm = vecnorm(processed_data.Magnetometer,2,2);
    plot(mag_norm)
%     ylim([-60 50])
    legend('L^{2}-norm')
    ylabel('\muT')
    title(sprintf('L^{2} norm (θ ≈ %d°)',atti(i)))

%     subplot(3,3,4+1*(i-1))
    subplot(3,4,3+4*(i-1))
    mag_x = processed_data.Magnetometer(:,1);
    mag_y = processed_data.Magnetometer(:,2);
    mag_z = processed_data.Magnetometer(:,3);
    plot([mag_x,mag_y,mag_z])
    ylabel('\muT')
%     hold on
%     plot(mag_x)    
%     plot(mag_y,'-.')    
%     plot(mag_z,':')    
%     hold off
%     grid on
    ylim([-80 40])
%     yticks(-80:40:50)
    lgd = legend('x','y','z','location','southeast');
    lgd.NumColumns = 3;
    title(sprintf('Raw 3-axis (θ ≈ %d°)',atti(i)))
    
%     subplot(3,3,7+1*(i-1))
    subplot(3,4,4+4*(i-1))
    rot_mag_x = rot_mag(:,1);
    rot_mag_y = rot_mag(:,2);
    rot_mag_z = rot_mag(:,3);
    plot([rot_mag_x,rot_mag_y,rot_mag_z])    
    ylabel('\muT')
%     hold on
%     plot(rot_mag_x)    
%     plot(rot_mag_y,'--')    
%     plot(rot_mag_z,':')    
%     hold off
%     grid on
    ylim([-80 40])   
%     yticks(-80:40:50)
    lgd = legend('x','y','z','location','southeast');
    lgd.NumColumns = 3;   
    title(sprintf('Calibrated 3-axis (θ ≈ %d°)',atti(i)))
end
%     ax = gca;
%     ax.LineStyleOrder = {'-','--',':'};
% set(gca,'LineStyleOrder',{'-','--o',':s'},'ColorOrder','g')
% sgtitle('Mag')

% sdf(gcf,'sj2')  % for response
sdf(gcf,'sj4')  % for paper
set(gcf,'units','points','position',[300,500,1600,700])
% tightfig
