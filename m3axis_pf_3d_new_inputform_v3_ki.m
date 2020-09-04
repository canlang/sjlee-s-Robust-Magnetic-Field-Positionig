clear; close all;
% purpose: convert the step length as a latent variable
% v3: add MaLoc 

% flag: Video save
video_flag = 0;

% success: 3(->path),4(<-path),6+,7+,8(loop-),9(loop+),10,11,12,13(-),14,16,17,18,19
% failure: 15,32,34,35,36
% loop: 43,

% t_input_idx = 63; ki
% t_input_idx = 29;
% t_input_idx = 36;

% switch t_input_idx
%     case 43
%         video_filename = 'm3axis_2d-space_pf_real-time-input_loop-traj';        %input index = 43;
%     otherwise
%         video_filename = 'm3axis_2d_pf_nonshiftedinput_newinput_rotated_6';%nonshifted input??
% end

% heading_noise = .01;        % heading noise candidates: .01, .50

method_name = {'ILoA', 'MaLoc'};
algo_idx = 2;
fprintf('%s testing...\n',method_name{algo_idx})

%%
% (1)
site_name = 'KI-1F';
map = magmap_construction('mats',site_name,.6);
lm.x = map(:,1);lm.y = map(:,2);
lM = map(:,3:5);

% (2)
% load('mats/magmap-ki-1f-0.6p.mat');
% lm.x = map(:,1);lm.y = map(:,2);
% lM = map(:,3:5);

% (3)
% data1 = readtable('sw_maps/dataset_ki_1f.csv');
% lM = [data1.mag_x,data1.mag_y,data1.mag_z];
% lm.x = data1.loc_x;lm.y = data1.loc_y;

% INTERPOLATION
% x = data1.loc_x;
% y = data1.loc_y;
% newlM = [];
% d = .7;
% for i=1:3
%     z = lM(:,i);
%     XI = min(x):d:max(x);
%     YI = min(y):d:max(y);
%     [X,Y] = meshgrid(XI,YI);
%     shp = alphaShape(x,y,5,'RegionThreshold',0);
%     in = inShape(shp,X,Y);
%     xg = X(in);
%     yg = Y(in);
%     zg = griddata(x,y,z,xg,yg,'nearest');         % 2. griddata() : INTERPOLATION
%     if isempty(newlM)
%         newlM = [xg,yg,zg];
%     else
%         newlM = [newlM,zg];
%         if ~isequal(newlM(:,1:2), [xg,yg])
%             disp 'error'
%         end
%     end
% end
% data1 = array2table(newlM(:,1:2), 'VariableNames',{'x','y'});
% map = [data1.x, data1.y, newlM(:,3:end)];
% lm.x = map(:,1);lm.y = map(:,2);
% lM = map(:,3:5);


% plot(lm.x,lm.y,'.')
%%
% #1. old testing data
% data2 = readtable('20171124 MagCoord3axisData.csv');
% #2. new collected data (realistic tracking)
tr_idx = 3;
% device_name = 'MATE20pro';
device_name = 'S9';

testfolder_name = sprintf('test-%s-%s',site_name,device_name);
% testfolder_name = 'test-KI-1F-S9';
test_path = sprintf('rawdata/%s',testfolder_name);
target_rawdata_paths = getNameFolds(test_path);
rawdata = load_rawdata(fullfile(test_path,target_rawdata_paths{tr_idx}));

% test_data_paths = dir('rawdata/test-ki-1f-iphone/*.csv');
% rawdata = load_rawdata(test_data_paths(tr_idx),'iPhone');
% 
% plot(rawdata.acc_norm)
% hold on
% plot(rawdata1.acc_norm)
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
site_name = 'KI-1F';
layout = jsondecode(fileread(sprintf('map/%s.json',site_name)));
% layout = loadjson(sprintf('map/%s.json',site_name));
% layout = loadjson('map/KI-1F.json');

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
% print -depsc2 env_setting.eps

% axis equal
% xlim([8 83])
% ylim([0 30])
% set(gcf,'units','points','position',[700,500,1500,700])

plot(shp,'FaceAlpha',.1,'EdgeColor','r')
legend('reference point','layout area')

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
eta = zeros(length(tM),n);
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
    
%     mu = [0,0];
% %     sigma = [0.10409786, 0.13461109; 0.13461109, 0.29744705];
% %     sigma = [1.43724175 -0.93884837;-0.93884837  0.74606608];
%     sigma = [0.02765426 0;0  0.04181993];
%     mvnRand = mvnrnd(mu,sigma,n);
%     ps.x = ps.x + cos(ps.mag_heading+euler(i,3)+1*pi/2).*ps.sl+ mvnRand(:,1);
%     ps.y = ps.y + sin(ps.mag_heading+euler(i,3)+1*pi/2).*ps.sl+ mvnRand(:,2);

    ps.x = ps.x + cos(ps.mag_heading+euler(i,3)).*ps.sl+ random('Uniform',-.1,.1,n,1);
    ps.y = ps.y + sin(ps.mag_heading+euler(i,3)).*ps.sl+ random('Uniform',-.1,.1,n,1);

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
    
% %     R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]/(rotMat(:,:,i))),0,...
% %         'UniformOutput',false);
% %     rotatedMag = cell2mat(cellfun(@(x)((x*tM(i,:)')'),R,'UniformOutput',false));
% 
% %     (1) UPDATE FUNCTION candidate 
% %     R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]*(rotMat(:,:,i))),ps.mag_heading,...
% %         'UniformOutput',false); 
% %     rotatedMag = cell2mat(cellfun(@(x)(x.'*tM(i,:)')',R,'UniformOutput',false));
% %     (2) UPDATE FUNCTION candidate  
% %     R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]/(rotMat(:,:,i))),ps.mag_heading,...
% %         'UniformOutput',false); 
% %     rotatedMag = cell2mat(cellfun(@(x)(x*tM(i,:)')',R,'UniformOutput',false));
% %     (3) UPDATE FUNCTION candidate  
%     R = arrayfun(@(x)(rotMat(:,:,i)*[cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1])...
%         ,-ps.mag_heading,'UniformOutput',false); 
%     rotatedMag = cell2mat(cellfun(@(x)(x.'*tM(i,:)')',R,'UniformOutput',false));
% 
% %     R = arrayfun(@(x) euler2rotMat(-euler(1,1),-euler(1,2),x),...
% %         -euler(1,3)-ps.mag_heading,'UniformOutput',false);
% %     rotatedMag = cell2mat(cellfun(@(x)(x*tM(i,:)')',R,'UniformOutput',false));
%     
%     % EUCLIDEAN
%     % TODO: may be more optimizable (DONE?maybe)
% %     mag_dist = diag(pdist2(rotatedMag,lM(I,:),'euclidean'));
%     mag_dist = sqrt(sum((rotatedMag-lM(I,:)).^2,2));
% %     mag_dist = bsxfun(@(x,y) pdist([x;y]), rotatedMag,lM(I,:));
%     
%     % COSINE
% %     mag_dist = diag(pdist2(rotatedMag,lM(I,:),'minkowski',3));

    if algo_idx == 1   
        % (1) tilting compensation
%         rotatedMag = getHeadingRotatedVector(ps.mag_heading, tM(i,:), rotMat(:,:,i));
        % (2) not applied tilting compenstation
        R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]),...
            ps.mag_heading+euler(i,3),'UniformOutput',false);
        rotatedMag = cell2mat(cellfun(@(x)((x*tM(i,:)')'),R,'UniformOutput',false));
        mag_dist = sqrt(sum((rotatedMag-lM(I,:)).^2,2));
    elseif algo_idx == 2
        observed = [vecnorm(lM(I,1:2),2,2), lM(I,3)];
        measured = [vecnorm(tM(i,1:2),2), tM(i,3)];
        mag_dist2 = pdist2(observed, measured,'mahalanobis',nearestSPD(cov(observed)));
        mag_dist = mag_dist2';
    else        % Unexpedted (wrong) case
        mag_dist = 0;
        disp('wrong feature algoithm input');
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
    eta(i,:) = ps.mag_heading';
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
    A = imread(sprintf('map/%s.png',site_name),'BackgroundColor',[1 1 1]);

%     scale = 0.05;
    xWorldLimits = [0 size(A,2)*scale];
    yWorldLimits = [0 size(A,1)*scale];
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
% figure
subplot(211)
% (1)
% histogram(wrapTo2Pi(ps.phy_heading))
% xlabel('Physical Heading')
% xlim([0 2*pi])
% set(gca,'XTick',0:pi/2:2*pi) 
% set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'}) 
% (2)
% plot(est(:,3))
plot(wrapTo2Pi(est(:,3)),'x-')
% (3)
% [lineh, bandsh] = fanChart(1:size(eta,1),eta, 'mean', 10:10:90, ...
%     'alpha', .2, 'colormap', {'shadesOfColor', [0 0 .8]});
% txt = strcat({'Pct'}, cellstr(int2str((20:20:80)')));
% legend([bandsh;lineh], [txt;{'Mean'}])

grid on
% ylim([-pi pi])
ylim([0 2*pi])
xlim([0,length(est)])
% set(gca,'YTick',-pi:pi/2:pi) 
% set(gca,'YTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi',}) 
set(gca,'YTick',0:pi/2:2*pi) 
set(gca,'YTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'}) 
xlabel('input index')
ylabel('\eta (rad)')

subplot(212)
% plot(euler(:,3))
% grid on
% ylim([-pi pi])
% xlim([0,length(est)])
% set(gca,'YTick',-pi:pi/2:pi) 
% set(gca,'YTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi',}) 

% histogram visualize
x = 1:size(eta,1);
y = eta;
xRep = repmat(x, 1, n);
h = histogram2(xRep(:),y(:),[size(eta,1), 100],'DisplayStyle','tile','Normalization','probability');
set(gca,'YTick',0:pi/2:2*pi) 
set(gca,'YTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'}) 

xlabel('input index')
ylabel('\eta (rad)')

% subplot(313)
% % stem(wrapTo2Pi(ps.mag_heading),ones(n,1))
% 
% % histogram(wrapTo2Pi(ps.mag_heading),n/2)
% % xlabel('Latest Time Magnetic Heading')
% % xlim([0 2*pi])
% % set(gca,'XTick',0:pi/2:2*pi) 
% % set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
% 
% turn_thr = .15;
% [~,turn_locs] = findpeaks(std_euler,'MinPeakHeight',turn_thr);
% findpeaks(std_euler,'MinPeakHeight',turn_thr)
% xlim([0,length(est)])


set(gcf,'units','points','position',[800,100,900,700])
sdf(gcf,'sj2')

% save('mats/iloa_success_ki_tr3.mat','est','eta')
%%
gt_filename = sprintf('est-result/ki-gt-s%d.json',tr_idx);
gt = jsondecode(fileread(gt_filename));
% converged_est = [est(converge_idx:end,1), est(converge_idx:end,2)];

% copied line 484
turn_thr = .15;
[~,turn_locs] = findpeaks(std_euler,'MinPeakHeight',turn_thr);

est_testpts = [est(turn_locs,1),est(turn_locs,2)];
disp('error to gt')
err2gt = diag(pdist2(gt(2:end,:), est_testpts));
disp(err2gt)
return
%% for papaer figure (without map image)
close all
err_std = std(err,0,2);
converge_idx = find(err_std <= 2,1);

axis equal

est2(:,1) = est(:,1)-12;
est2(:,2) = est(:,2)-21;


hold on
plot(est2(1:converge_idx,1), est2(1:converge_idx,2),'-xr')
plot(est2(converge_idx:end,1), est2(converge_idx:end,2),'-b')
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
% xlim([0.6900   34.1538]);
% ylim([-15.0702   10.5512]);

set(gcf,'units','points','position',[500,500,800,600])
sdf(gcf,'sj2')
% print -depsc2 eps/18times_repeat_circle_raw.eps
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