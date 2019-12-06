clear; close all;
% purpose: convert the step length as a latent variable

% flag: Video save
video_flag = 0;

% success: 3(->path),4(<-path),6+,7+,8(loop-),9(loop+),10,11,12,13(-),14,16,17,18,19
% failure: 15,32,34,35,36
% loop: 43,

t_input_idx = 3;
% t_input_idx = 29;
% t_input_idx = 36;

switch t_input_idx
    case 43
        video_filename = 'm3axis_2d-space_pf_real-time-input_loop-traj';        %input index = 43;
    otherwise
        video_filename = 'm3axis_2d_pf_nonshiftedinput_newinput_rotated_6';%nonshifted input??
end

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
[pks,locs] = findpeaks(processed_data.acc_norm,'MinPeakDistance',...
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

A = imread('N1-7F.png','BackgroundColor',[1 1 1]);

xWorldLimits = [-1 1650/20];
yWorldLimits = [-1 660/20];
RA = imref2d(size(A),xWorldLimits,yWorldLimits);
imshow(flipud(A),RA);
axis xy;


% draw learning data
hold on
plot(lm.x,lm.y,'.','MarkerSize', 10)
% plot(lM(:,1),lM(:,2),'.','MarkerSize', 10)
% for save eps
legend('reference point')
sdf(gcf,'sj2')
% print -depsc2 env_setting.eps

% axis equal
% xlim([8 83])
% ylim([0 30])
% set(gcf,'units','points','position',[700,500,1500,700])
plot(shp,'FaceAlpha',.5,'EdgeColor','r')


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

% draw particle
hold on 
% h_ps = plot(ps.x,ps.y,'.','MarkerSize',8);
h_ps = scatter(ps.x,ps.y,20,'c','filled','MarkerFaceAlpha',.2);
h_pm = plot(mean(ps.x),mean(ps.y),'ms');
% h_gt = plot(data2.coord_x(1),data2.coord_y(1),'s','MarkerSize',10,...
%     'MarkerEdgeColor','b',...
%     'MarkerFaceColor',[0.5,0.5,0.5]);
hold off


% test data matrix
% tM = [data2.mag_x,data2.mag_y,data2.mag_z];
% rotate test data 
% x = -pi/2;
x = 0;
R = [cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1];

% tM = (R*tM')';
% tM = (R*tM')'-25;
% 25 offset meaning TestData shift: because magnetometer's calibration not matched
tM = processed_data.Magnetometer(locs,:);

% test result matrix
est = zeros(length(tM),3);
err = zeros(length(tM),n);
%%
switch video_flag
    case 1
        v = VideoWriter(strcat('movie_files/',video_filename,'.mp4'),'MPEG-4');
        v.FrameRate = 10;
        v.Quality = 100;
        open(v);
        frame = getframe(gcf);
        writeVideo(v,frame);
end

for i = 1:length(tM)
    % ================ PREDICTION
%     ps.x = bsxfun(@(x,y) x + cos(y),ps.x,ps.phy_heading);
%     ps.y = bsxfun(@(x,y) x + sin(y),ps.y,ps.phy_heading);
%     ps.x = ps.x + cos(ps.mag_heading-pi/2);     % when heading shifted input date
%     ps.y = ps.y + sin(ps.mag_heading-pi/2);
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
    
    mu = [0,0];
%     sigma = [0.10409786, 0.13461109; 0.13461109, 0.29744705];
    sigma = [2.97925885  0.02243997;0.02243997  0.98844656];
    mvnRand = mvnrnd(mu,sigma,n);
    ps.x = ps.x + cos(ps.mag_heading+euler(i,3)).*ps.sl+ mvnRand(:,1);
    ps.y = ps.y + sin(ps.mag_heading+euler(i,3)).*ps.sl+ mvnRand(:,2);

%     ps.x = ps.x + cos(ps.mag_heading+euler(i,3)).*ps.sl+ random('Uniform',-.1,.1,n,1);
%     ps.y = ps.y + sin(ps.mag_heading+euler(i,3)).*ps.sl+ random('Uniform',-.1,.1,n,1);

%     ps.x = ps.x + cos(ps.mag_heading+euler(i,3))*sl;
%     ps.y = ps.y + sin(ps.mag_heading+euler(i,3))*sl;
%     ps.x = ps.x + cos(ps.mag_heading+euler(i,3))*sl + random('Uniform',-1,1,n,1);
%     ps.y = ps.y + sin(ps.mag_heading+euler(i,3))*sl + random('Uniform',-1,1,n,1);
%     ps.x = ps.x + cos(ps.mag_heading+d_heading(i))*sl + random('Uniform',-1,1,n,1);
%     ps.y = ps.y + sin(ps.mag_heading+d_heading(i))*sl + random('Uniform',-1,1,n,1);
%     ps.x = ps.x + ps.stlng.*cos(ps.heading);
%     ps.y = ps.y + ps.stlng.*sin(ps.heading);
    
    % ================ UPDATE    
    % 1. find (geo-locational) nearest learning data
    [phy_dist,I] = pdist2([lm.x,lm.y],[ps.x,ps.y],'euclidean','Smallest',1);
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
    R = arrayfun(@(x)(rotMat(:,:,i)*[cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]),ps.mag_heading,...
        'UniformOutput',false); 
    rotatedMag = cell2mat(cellfun(@(x)(x.'*tM(i,:)')',R,'UniformOutput',false));
    
    % EUCLIDEAN
    % TODO: may be more optimizable (DONE?maybe)
%     mag_dist = diag(pdist2(rotatedMag,lM(I,:),'euclidean'));
    mag_dist = sqrt(sum((rotatedMag-lM(I,:)).^2,2));
%     mag_dist = bsxfun(@(x,y) pdist([x;y]), rotatedMag,lM(I,:));
    
    % COSINE
%     mag_dist = diag(pdist2(rotatedMag,lM(I,:),'minkowski',3));


    if all(mag_dist) == 0
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
%     phy_move_noise_range = 2;
%     ps.x = ps.x(resample_idx) + phy_move_noise_range*rand(n,1) - phy_move_noise_range/2;
%     ps.y = ps.y(resample_idx) + phy_move_noise_range*rand(n,1) - phy_move_noise_range/2;
%     ps.x = ps.x(resample_idx) + random('normal',0,.5,n,1);
%     ps.y = ps.y(resample_idx) + random('normal',0,.5,n,1);
    ps.x = ps.x(resample_idx);
    ps.y = ps.y(resample_idx);
    ps.mag_heading = ps.mag_heading(resample_idx);
    ps.sl = ps.sl(resample_idx);
    
%     if std_euler(i) > 0.03
%         ps.mag_heading = ps.mag_heading(resample_idx)+random('normal',0,std_euler(i),n,1);
%     else
%         ps.mag_heading = ps.mag_heading(resample_idx);
%     end
%     ps.mag_heading = ps.mag_heading(resample_idx)+random('normal',0,heading_noise,n,1);    
%     ps.mag_heading = ps.mag_heading(resample_idx);
%     ps.phy_heading = ps.phy_heading(resample_idx)+random('normal',0,.001,n,1);

%     ps.heading = ps.heading(resample_idx)+random('Uniform', -pi/10,pi/10,n,1);
%     ps.stlng = ps.stlng(resample_idx) + random('normal',0,1,n,1);
%     ps.prob = init_prob;
   
    set(h_ps,'XData',ps.x,'YData',ps.y)                             % ps result
    set(h_pm,'XData',mean(ps.x),'YData',mean(ps.y));
%     set(h_gt,'XData',data2.coord_x(i),'YData',data2.coord_y(i))     % ground truth
    drawnow
    
%     est(i,:) = [mean(ps.x),mean(ps.y),mean(angdiff(ps.mag_heading,0))];
%     est(i,:) = [mean(ps.x),mean(ps.y),mean(abs(angdiff(ps.mag_heading,0)))];        %x,y,heading
    est(i,:) = [mean(ps.x),mean(ps.y),circ_mean(ps.mag_heading)]; 
    
%     err(i,:) = pdist2([data2.coord_x(i),data2.coord_y(i)],[ps.x,ps.y]);
    err(i,:) = pdist2([mean(ps.x),mean(ps.y)],[ps.x,ps.y]);
%     break
%     pause(.1)
    switch video_flag
        case 1
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end
switch video_flag, case 1, close(v);end
% %%
% close all
% subplot(121)
% err_std = std(err,0,2);
% boxplot(find(err_std < 2,1))
% subplot(122)
% est_err = sqrt(sum((est(:,1:2)-[data2.coord_x data2.coord_y]).^2,2));
% cdfplot(est_err(find(err_std < 2,1):end))
% 
% set(gcf,'units','points','position',[200,500,2000,800])
% sdf(gcf,'sj2')
% return

%%
figure
err_std = std(err,0,2);
converge_idx = find(err_std <= 2,1);
tracking_show_with_map = true;
if tracking_show_with_map
    A = imread('N1-7F.png','BackgroundColor',[1 1 1]);

    xWorldLimits = [-1 1650/20];
    yWorldLimits = [-1 660/20];
    RA = imref2d(size(A),xWorldLimits,yWorldLimits);
    imshow(flipud(A),RA);
    axis xy;

    axis equal
end

hold on
plot(est(1:converge_idx,1), est(1:converge_idx,2),'xr')
plot(est(converge_idx:end,1), est(converge_idx:end,2),'+b')
legend('\sigma>2','\sigma\leq2')
% fprintf('Mean Err.=%.2f\n',est(converge_idx:end,1), est(converge_idx:end,2));
% cd = [uint8(jet(n)*255) uint8(ones(n,1))].';
% 
% drawnow
% set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cd)

sdf(gcf,'sj')
% print -depsc2 eps/18times_repeat_circle.eps

%%
figure
subplot(311)
% (1)
% histogram(wrapTo2Pi(ps.phy_heading))
% xlabel('Physical Heading')
% xlim([0 2*pi])
% set(gca,'XTick',0:pi/2:2*pi) 
% set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'}) 
% (2)
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
% stem(wrapTo2Pi(ps.mag_heading),ones(n,1))
histogram(wrapTo2Pi(ps.mag_heading),n/2)
xlabel('Latest Time Magnetic Heading')
xlim([0 2*pi])
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})

set(gcf,'units','points','position',[800,100,900,700])
sdf(gcf,'sj')

return
%% for papaer figure (without map image)
close all
err_std = std(err,0,2);
converge_idx = find(err_std <= 2,1);

axis equal

est2(:,1) = est(:,1)-12;
est2(:,2) = est(:,2)-21;


hold on
plot(est2(1:converge_idx,1), est2(1:converge_idx,2),'xr')
plot(est2(converge_idx:end,1), est2(converge_idx:end,2),'+b')
legend('\sigma>2','\sigma\leq2')
hold off
% set(gca,'xticklabel',[])
% set(gca,'yticklabel',[])
% xticks('auto')
% yticks('auto')
pbaspect([1 1 1])
grid on
grid minor
% set(gca,'xtick',[])
% set(gca,'ytick',[])
xlabel('x (m)');
ylabel('y (m)');
xlim([0.6900   34.1538]);
ylim([-15.0702   10.5512]);

set(gcf,'units','points','position',[500,500,800,600])
sdf(gcf,'sj2')
print -depsc2 eps/18times_repeat_circle_raw.eps
return



%% positioning and heading error (x: location, y: err)
figure
est_err = sqrt(sum((est(:,1:2)-[data2.coord_x data2.coord_y]).^2,2));
subplot(211)
plot(est_err,'s-')
xlabel('Index of Test Data ')
ylabel('Error Distance (m)')
subplot(212)
plot(abs(angdiff(wrapTo2Pi(est(:,3)),0)))
% ylim([0 2*pi])
xlabel('Index of Test Data ')
ylabel('Heading Error (rad)')
set(gcf,'units','points','position',[300,100,1000,800])
sdf(gcf,'sj2')
%%
clf
err = sqrt(sum((est(:,1:2)-[data2.coord_x data2.coord_y]).^2,2));
yyaxis left
plot(err,'s:','MarkerSize',8)
xlabel('Index of Test Data ')
ylabel('Error Distance (m)')

yyaxis right
plot(abs(angdiff(wrapTo2Pi(est(:,3)),0)),'*-.','MarkerSize',8)
% ylim([0 2*pi])
ylabel('Heading Error (rad)')
legend('Error of Position','Error of Heading')
set(gcf,'units','points','position',[1300,500,800,300])
sdf(gcf,'sj2')


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