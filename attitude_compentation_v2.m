clearvars;close all;clc;
target_rawdata_paths = getNameFolds('rawdata');
% data_idxes = 87:89;
% data_idxes = 92:-1:90;      % 0 45 84
% data_idxes = 92:94;      % 0 20 90
data_idxes = [92,91,94];     % 0 45 90
% data_idxes = 95:-1:93;
atti = [0,45,90];

mag_norm = cell(1,3);
raw_mags = cell(1,3);
rot_mags = cell(1,3);

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
    mag_norm{i} = vecnorm(processed_data.Magnetometer,2,2);
    raw_mags{i} = processed_data.Magnetometer;
    rot_mags{i} = rot_mag;
end
%%
close all; 
n_color = 3;
bar_cmap = cool(n_color);
set(0,'DefaultAxesColorOrder',bar_cmap);

ax_LineStyleOrder = {'-','--',':'};
for i=1:3
    subplot(3,4,1+4*(i-1))
    A = imread(sprintf('jpg/attitude_%d.jpg',atti(i)));    
    R = imrotate(A,-90);
    I = imcrop(R,[0 600 3024 2800]);
    imshow(I)
    title(sprintf('Atti.Group#%d (θ ≈ %d°)',i,atti(i)))
    
    subplot(3,4,6)    
    hold on
    plot(mag_norm{i},ax_LineStyleOrder{i})
%     ylim([-60 50])
    legend('Atti.G.#1','Atti.G.#2','Atti.G.#3');
    ylabel('\muT')
    title(sprintf('L^{2}-norm',atti(i)))
    
    ylabels = {'x','y','z'};
    
    for j=1:3
        subplot(3,4,3+(j-1)*4)
        hold on
        raw_mag = raw_mags{i}(:,j);
        plot(raw_mag,ax_LineStyleOrder{i})
        ylabel('\muT')
        title(sprintf('raw %s',ylabels{j}));
%         legend('Atti.G.#1','Atti.G.#2','Atti.G.#3');
    end
    
    for j=1:3
        subplot(3,4,4+(j-1)*4)
        hold on
        rot_mag = rot_mags{i}(:,j);
        plot(rot_mag,ax_LineStyleOrder{i})
        ylabel('\muT')
        title(sprintf('calibrated %s',ylabels{j}));
%         legend('Atti.G.#1','Atti.G.#2','Atti.G.#3');
    end
end
set(gcf,'units','points','position',[300,500,1100,500])

%     ax = gca;
%     ax.LineStyleOrder = {'-','--',':'};
% set(gca,'LineStyleOrder',{'-','--o',':s'},'ColorOrder','g')
% sgtitle('Mag')

% sdf(gcf,'sj2')  % for response
% sdf(gcf,'sj5')  % for paper
sdf(gcf,'paperwidth')  % for paper

% tightfig
