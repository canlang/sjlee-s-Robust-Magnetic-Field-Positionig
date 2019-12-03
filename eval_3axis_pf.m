function traj_errors = eval_3axis_pf(a,b,c)
% clear; close all;
% purpose: evaluate score the trajectory
% from: revised from
% m3axis_particle_filter_3d_std_validation_newinput_form_v2.m file

t_input_idx = 3;
% t_input_idx = 29;
% t_input_idx = 36;

% heading_noise = .01;        % heading noise candidates: .01, .50
%%
map = magmap_construction('mats',.5);
lm.x = map(:,1);lm.y = map(:,2);
lM = map(:,3:5);

%%
% #1. old testing data
% data2 = readtable('20171124 MagCoord3axisData.csv');
% #2. new collected data (realistic tracking)
target_rawdata_paths = getNameFolds('rawdata');
rawdata = load_rawdata(fullfile('rawdata',target_rawdata_paths{t_input_idx}));

%% resample for synchronize
rate = 2e-2;
processed_data = resample_rawdata(rawdata,rate);

%% find step point (step) and time labeling
% threshold should be tuned experimentally to match a person's level 
minPeakHeight = std(processed_data.acc_norm);       
[~,locs] = findpeaks(processed_data.acc_norm,'MinPeakDistance',...
    .3/rate,'MinPeakHeight',minPeakHeight);   % .3s 이내의 피크는 무시 (가정: 발걸음이 .3초 내에는 2번이뤄지지 않음)\

%%
addpath(genpath('madgwick_algorithm_matlab'));
AHRS = MadgwickAHRS('SamplePeriod', rate, 'Beta', 0.1); % sample rate: 2e-2

quaternion = zeros(length(processed_data.Time), 4);
for t = 1:length(processed_data.Time)
    AHRS.Update(processed_data.Gyroscope(t,:) * (pi/180), ...
        processed_data.Accelerometer(t,:), ...
        processed_data.Accelerometer(t,:));	% gyroscope units must be radians
    quaternion(t, :) = AHRS.Quaternion;
end

% euler : use conjugate for sensor frame relative to Earth and convert to degrees.
euler = quatern2euler(quaternConj(quaternion(locs,:)));	
% euler = quatern2euler(quaternConj(quaternion)) * (180/pi);
d_heading = [0;diff(euler(:,3))];

rotMat = quatern2rotMat(quaternion(locs,:));
% rotMat = quatern2rotMat(quaternConj(quaternion(locs,:)));

std_euler = stdfilt(unwrap(euler(:,3)));
% std_euler = stdfilt(deg2rad(unwrap(euler(:,3))));
% plot(processed_data.Time, std_euler)

%% inbound & outbound
layout = loadjson('N1-7F-HiRes2.json');

x = layout.in(:,1);
y = layout.in(:,2);
shp = polyshape(x,y);
for i = 1:length(layout.out)
    ox = layout.out{i}(:,1);
    oy = layout.out{i}(:,2);
    shp = subtract(shp,polyshape(ox,oy));
end

%%
% initialize particle
n = 2000;
% 1. only road
rand_idx = randi(length(lm.x),n,1);
ps.x = lm.x(rand_idx);
ps.y = lm.y(rand_idx);

% 2. all area
% ps.x = random('Uniform', min(lm.x),max(lm.x),n,1);
% ps.y = random('Uniform', min(lm.y),max(lm.y),n,1);

% 3. initial area
% ps.x = data2.coord_x(1)+random('normal',0,5,n,1);
% ps.y = data2.coord_y(1)+random('normal',0,5,n,1);
% ps.y = 15+random('normal',0,5,n,1);

ps.sl = ones(n,1)*.7;
ps.mag_heading = random('Uniform', 0,2*pi,n,1);
% ps.phy_heading = random('Uniform', 0,2*pi,n,1);
ps.prob = ones(n,1)*(1/n);
% ps.stlng = ones(n,1) + random('Uniform', -.1,.1,n,1);

x = 0;
R = [cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1];

tM = processed_data.Magnetometer(locs,:);

% test result matrix
est = zeros(length(tM),3);
err = zeros(length(tM),n);
%%

for i = 1:length(tM)
    % ================ PREDICTION    
    if i>1
%         Halpha = pi-est(i-1,3);
        Halpha = 0;
        ps.mag_heading = ps.mag_heading+random('normal',Halpha,.03,n,1);     % candidate, .08
    else
%     ps.mag_heading = ps.mag_heading+euler(i,3);
        ps.mag_heading = ps.mag_heading+random('normal',0,.03,n,1);     % candidate, .08
    end
    
    ps.sl = .7 + random('normal',0,.5,n,1);
    
    mu = [0,0];
    sigma = [a b; b c];
    mvnRand = mvnrnd(mu,sigma,n);
    ps.x = ps.x + cos(ps.mag_heading+euler(i,3)).*ps.sl+ mvnRand(:,1);
    ps.y = ps.y + sin(ps.mag_heading+euler(i,3)).*ps.sl+ mvnRand(:,2);
%     ps.x = ps.x + cos(ps.mag_heading+euler(i,3)).*ps.sl+ random('Uniform',-.1,.1,n,1);
%     ps.y = ps.y + sin(ps.mag_heading+euler(i,3)).*ps.sl+ random('Uniform',-.1,.1,n,1);
    
    % ================ UPDATE    
    % 1. find (geo-locational) nearest learning data
    [phy_dist,I] = pdist2([lm.x,lm.y],[ps.x,ps.y],'euclidean','Smallest',1);
    % 2. calculate Rotated magnetic field data and magnetic distance 
    R = arrayfun(@(x)(rotMat(:,:,i)*[cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]),ps.mag_heading,...
        'UniformOutput',false); 
    rotatedMag = cell2mat(cellfun(@(x)(x.'*tM(i,:)')',R,'UniformOutput',false));
    
    % EUCLIDEAN
    mag_dist = sqrt(sum((rotatedMag-lM(I,:)).^2,2));
    
    % COSINE
%     mag_dist = diag(pdist2(rotatedMag,lM(I,:),'minkowski',3));

    if ~all(mag_dist)
        break
    end
    ps.prob = 1./(mag_dist);
    in = isinterior(shp,ps.x,ps.y);
%     in = inShape(shp,[ps.x,ps.y]);
    ps.prob(~in) = 0;
    ps.prob(phy_dist'>3) = 0;
    
    if sum(ps.prob) == 0
        rand_idx = randi(length(lm.x),n,1);
        ps.x = lm.x(rand_idx);
        ps.y = lm.y(rand_idx);
        ps.mag_heading = random('Uniform', 0,2*pi,n,1);
%         ps.phy_heading = random('Uniform', 0,2*pi,n,1);
        ps.prob = ones(n,1)*(1/n);
    else
        ps.prob = ps.prob./sum(ps.prob);
    end
        
    % ================ RESAMPLE
    resample_idx = randsample(1:n,n,true,ps.prob);
    ps.x = ps.x(resample_idx);
    ps.y = ps.y(resample_idx);
    ps.mag_heading = ps.mag_heading(resample_idx);
    ps.sl = ps.sl(resample_idx);

    est(i,:) = [mean(ps.x),mean(ps.y),circ_mean(ps.mag_heading)];     
end

%%

err_std = std(err,0,2);
converge_idx = find(err_std <= 2,1);
% [9,21],[71,21]
gt_x = linspace(9,70.7,length(tM));
gt_y = 21.2*ones(1,length(tM));

traj_errors = diag(pdist2([est(converge_idx:end,1),est(converge_idx:end,2)],...
    [gt_x(converge_idx:end)',gt_y(converge_idx:end)']));

    %% local functions -----------------------------------------------------------------
    % 1st local function
    function nameFolds = getNameFolds(pathFolder)
        d = dir(pathFolder);
        isub = [d(:).isdir]; %# returns logical vector
        nameFolds = {d(isub).name}';
        nameFolds(ismember(nameFolds,{'.','..'})) = [];
        % nameFolds(~cellfun(@isempty,regexp(nameFolds, '\d{6}'))) = [];
    end

    % 2nd local function
    function data = resample_rawdata(rawdata,rate)
        % Var1 : 3 vectors of Accelerometer, Var2: norm vector of Acc.
        T_acc = timetable(seconds(rawdata.acc(:,2)/1e9),rawdata.acc(:,3:5),rawdata.acc_norm,...
            'VariableNames',{'acc','acc_norm'});
        % Var1 : 3 vectors of gyroscope
        T_gyr = timetable(seconds(rawdata.gyr(:,2)/1e9),rawdata.gyr(:,3:5),...
            'VariableNames',{'gyr'});
        T_mag = timetable(seconds(rawdata.mag(:,2)/1e9),rawdata.mag(:,3:5),...
            'VariableNames',{'mag'});
        T_acc = sortrows(T_acc);
        T_gyr = sortrows(T_gyr);
        T_mag = sortrows(T_mag);

        TT = synchronize(T_acc,T_gyr,T_mag,'regular','linear','TimeStep',seconds(rate));

        data.Accelerometer = TT.acc;
        data.Gyroscope = TT.gyr*(180/pi);
        data.Magnetometer = TT.mag;
        data.acc_norm = TT.acc_norm;

        data.Time = seconds(TT.Time(:));
        % data.Time = seconds(TT.Time(:)-(TT.Time(1)));
        data.Rate = median(diff(data.Time)); % cal sample rate
    end
end
