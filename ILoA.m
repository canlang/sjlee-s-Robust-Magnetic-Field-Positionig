function err = ILoA(site_name,device_name,tr_idx,intp_intv,vis_flag,dist_idx, varargin)
% purpose: convert the localization function - for evaluation loop

%%
if nargin < 6
    dist_idx = 1;
end

map = loadMagneticMap('mats',site_name,intp_intv);
lm.x = map(:,1);lm.y = map(:,2);
lM = map(:,3:5);

testfolder_name = sprintf('test-%s-%s',site_name,device_name);
test_path = sprintf('rawdata/%s',testfolder_name);
if isequal(device_name, 'iphone')
    test_rawdata_paths = getNameFiles(test_path);
    rawdata = load_rawdata(fullfile(test_path,test_rawdata_paths{tr_idx}),'iPhone');    % iOS
else        
    test_rawdata_paths = getNameFolds(test_path);
    rawdata = load_rawdata(fullfile(test_path,test_rawdata_paths{tr_idx}));             % Android
end


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

rotMat = quatern2rotMat(quaternion(locs,:));
% rotMat = quatern2rotMat(quaternConj(quaternion(locs,:)));

std_euler = stdfilt(unwrap(euler(:,3)));
% [~,turn_locs] = findpeaks(std_euler,'MinPeakHeight',.3);              %
% for noise optimization, during turning area
%% inbound & outbound
% site_name = 'KI-1F';
layout = jsondecode(fileread(sprintf('map/%s.json',site_name)));

x = layout.in(:,1);
y = layout.in(:,2);
shp = polyshape(x,y);
if ~isempty(layout.out)
    for i = 1:length(layout.out)
        ox = layout.out{i}(:,1);
        oy = layout.out{i}(:,2);
        shp = subtract(shp,polyshape(ox,oy));
    end
end

if vis_flag
    A = imread(sprintf('map/%s.png',site_name),'BackgroundColor',[1 1 1]);
    scale = 0.05;
    xWorldLimits = [0 size(A,2)*scale];
    yWorldLimits = [0 size(A,1)*scale];
    % xWorldLimits = [-1 1650*.05];
    % yWorldLimits = [-1 660*.05];
    RA = imref2d(size(A), xWorldLimits, yWorldLimits);
    imshow(flipud(A),RA);
    axis xy;


    % draw learning data
    hold on
    plot(lm.x,lm.y,'.','MarkerSize', 10)
    % plot(lM(:,1),lM(:,2),'.','MarkerSize', 10)
    % for save eps
    % legend('reference point')
    sdf(gcf,'sj2')

    plot(shp,'FaceAlpha',.1,'EdgeColor','r')
    legend('reference point','layout area')
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
% ps.x = 45+random('normal',0,1,n,1);
% ps.y = 45+random('normal',0,1,n,1);

ps.sl = ones(n,1)*.7;
ps.mag_heading = random('Uniform', 0,2*pi,n,1);
ps.prob = ones(n,1)*(1/n);


if vis_flag
    hold on 
    % h_ps = plot(ps.x,ps.y,'.','MarkerSize',8);
    h_ps = scatter(ps.x,ps.y,20,'c','filled','MarkerFaceAlpha',.2);
    h_pm = plot(mean(ps.x),mean(ps.y),'ms');
    % h_gt = plot(data2.coord_x(1),data2.coord_y(1),'s','MarkerSize',10,...
    %     'MarkerEdgeColor','b',...
    %     'MarkerFaceColor',[0.5,0.5,0.5]);
    hold off
end

tM = processed_data.Magnetometer(locs,:);

est = zeros(length(tM),3);
%%

for i = 1:length(tM)
    % ================ PREDICTION 
    hN = .03;      % heading noise value
    if i>1
%         Halpha = pi-est(i-1,3);
        Halpha = 0;
        ps.mag_heading = ps.mag_heading+random('normal',Halpha,hN,n,1);     % candidate, .08
    else
%     ps.mag_heading = ps.mag_heading+euler(i,3);
        ps.mag_heading = ps.mag_heading+random('normal',0,hN,n,1);     % candidate, .08
    end
    
    ps.sl = .7 + random('normal',0,.5,n,1);

    ps.x = ps.x + cos(ps.mag_heading+euler(i,3)).*ps.sl+ random('Uniform',-.1,.1,n,1);
    ps.y = ps.y + sin(ps.mag_heading+euler(i,3)).*ps.sl+ random('Uniform',-.1,.1,n,1);
    
    % ================ UPDATE    
    % 1. find (geo-locational) nearest learning data
    [phy_dist,I] = findNearestLocation([lm.x,lm.y],[ps.x,ps.y]);
    % 2. calculate Rotated magnetic field data and magnetic distance    
%     R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]/(rotMat(:,:,i))),0,...
%         'UniformOutput',false);
%     rotatedMag = cell2mat(cellfun(@(x)((x*tM(i,:)')'),R,'UniformOutput',false));

%     (1) UPDATE FUNCTION candidate 
%     R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]*(rotMat(:,:,i))),ps.mag_heading,...
%         'UniformOutput',false); 
%     rotatedMag = cell2mat(cellfun(@(x)(x.'*tM(i,:)')',R,'UniformOutput',false));
%     (2) UPDATE FUNCTION candidate  
%     R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]/(rotMat(:,:,i))),ps.mag_heading,...
%         'UniformOutput',false); 
%     rotatedMag = cell2mat(cellfun(@(x)(x*tM(i,:)')',R,'UniformOutput',false));
%     (3) UPDATE FUNCTION candidate  
    if dist_idx == 1
        % TODO: have to check other rotated version. 
%         R = arrayfun(@(x)(rotMat(:,:,i)*[cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1])...
%             ,-ps.mag_heading,'UniformOutput',false); 
%         rotatedMag = cell2mat(cellfun(@(x)(x.'*tM(i,:)')',R,'UniformOutput',false));
        rotatedMag = getHeadingRotatedVector(ps.mag_heading, tM(i,:), rotMat(:,:,i));

    %     R = arrayfun(@(x) euler2rotMat(-euler(1,1),-euler(1,2),x),...
    %         -euler(1,3)-ps.mag_heading,'UniformOutput',false);
    %     rotatedMag = cell2mat(cellfun(@(x)(x*tM(i,:)')',R,'UniformOutput',false));

        % EUCLIDEAN
        % TODO: may be more optimizable (DONE?maybe)
    %     mag_dist = diag(pdist2(rotatedMag,lM(I,:),'euclidean'));
        mag_dist = sqrt(sum((rotatedMag-lM(I,:)).^2,2));
    %     mag_dist = bsxfun(@(x,y) pdist([x;y]), rotatedMag,lM(I,:));

        % COSINE
    %     mag_dist = diag(pdist2(rotatedMag,lM(I,:),'minkowski',3));
    elseif dist_idx == 2   % related work: Horizontal & Vertical components, MaLoc, B_h,B_v    
        observed = [vecnorm(lM(I,1:2),2,2), lM(I,3)];
        measured = [vecnorm(tM(i,1:2),2), tM(i,3)];
        mag_dist2 = pdist2(observed, measured,'mahalanobis',nearestSPD(nancov(observed)));
        mag_dist = mag_dist2';
    else
        disp('ILoA.m - wrong method input arg')
        mag_dist = 0;
    end

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
    if vis_flag
        set(h_ps,'XData',ps.x,'YData',ps.y)                             % ps result
        set(h_pm,'XData',mean(ps.x),'YData',mean(ps.y));
        drawnow        
    end
end

turn_thr = .1;
[~,turn_locs] = findpeaks(std_euler,'MinPeakHeight',turn_thr,'MinPeakDistance',5);
    
if vis_flag
    figure
    subplot(311)
    plot(est(:,3))
    grid on
    ylim([-pi pi])
    set(gca,'YTick',-pi:pi/2:pi) 
    set(gca,'YTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi',}) 
    xlabel('input index')
    ylabel('Magnetic \Delta (rad)')

    subplot(312)
    plot(euler(:,3))
    grid on
    ylim([-pi pi])
    set(gca,'YTick',-pi:pi/2:pi) 
    set(gca,'YTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi',}) 
    xlabel('input index')
    ylabel('Relative \psi (rad)')

    subplot(313)
    % % stem(wrapTo2Pi(ps.mag_heading),ones(n,1))
    % histogram(wrapTo2Pi(ps.mag_heading),n/2)
    % xlabel('Latest Time Magnetic Heading')
    % xlim([0 2*pi])
    % set(gca,'XTick',0:pi/2:2*pi) 
    % set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
    findpeaks(std_euler,'MinPeakHeight',turn_thr,'MinPeakDistance',5)

    set(gcf,'units','points','position',[800,100,900,700])
    sdf(gcf,'sj')
end

%%
gt_filename = sprintf('est-result/ki-gt-s%d.json',tr_idx);
gt = jsondecode(fileread(gt_filename));
% converged_est = [est(converge_idx:end,1), est(converge_idx:end,2)];
est_testpts = [est(turn_locs,1),est(turn_locs,2)];
err = diag(pdist2(gt(2:end,:), est_testpts));
end
%% local functions -----------------------------------------------------------------
% 1. local function
function data = resample_rawdata(rawdata,rate)
% Var1 : 3 vectors of Accelerometer, Var2: norm vector of Acc.
% seconds(rawdata.acc(:,2)/1e9)
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

% TT = synchronize(T_acc,T_gyr,T_mag,'regular','nearest','TimeStep',seconds(rate));
TT = synchronize(T_acc,T_gyr,T_mag,'regular','linear','TimeStep',seconds(rate));

data.Accelerometer = TT.acc;
data.Gyroscope = TT.gyr*(180/pi);
data.Magnetometer = TT.mag;
data.acc_norm = TT.acc_norm;

data.Time = seconds(TT.Time(:));
% data.Time = seconds(TT.Time(:)-(TT.Time(1)));
data.Rate = median(diff(data.Time)); % cal sample rate
end