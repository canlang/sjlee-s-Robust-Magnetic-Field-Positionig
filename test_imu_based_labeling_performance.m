clearvars;close all;
% most source referenced from ILoA.m
% load('mats/magmap-n1-7f-stepw-full-opt.mat','map');
load('mats/magmap-n1-7f-stepw-t-turn.mat','map');

lm.x = map(:,1);lm.y = map(:,2);
lM = map(:,3:5);

test_dataset = '200715_160602_093_n1_r_turn_return';
% test_dataset = '190228_161229_230-B3C';
% test_dataset = '190228_160953_233-D431A';
% test_dataset = '190228_shoon_A13B';


test_path = sprintf('rawdata/%s',test_dataset);
rawdata = load_rawdata(test_path);

v = VideoWriter(sprintf('vids/labeling-test/%s.mp4',test_dataset),'MPEG-4');
open(v);

%% resample for synchronize
rate = 2e-2;
testdata = resample_rawdata2(rawdata,rate);

minPeakHeight = std(testdata.acc_norm);       
[~,locs] = findpeaks(testdata.acc_norm,'MinPeakDistance',.3/rate,'MinPeakHeight',minPeakHeight);

addpath(genpath('madgwick_algorithm_matlab'));
AHRS = MadgwickAHRS('SamplePeriod', rate, 'Beta', 0.1); % sample rate: 2e-2

quaternion = zeros(length(testdata.Time), 4);
for t = 1:length(testdata.Time)
    AHRS.Update(testdata.Gyroscope(t,:),testdata.Accelerometer(t,:),testdata.Accelerometer(t,:));
    quaternion(t, :) = AHRS.Quaternion;
end
euler = quatern2euler(quaternConj(quaternion(locs,:)));	
rotMat = quatern2rotMat(quaternion(locs,:));

%% particle filter
algorithm = 2;
n = 3000;

rand_idx = randi(length(lm.x),n,1);
ps.x = lm.x(rand_idx);
ps.y = lm.y(rand_idx);
ps.sl = ones(n,1)*.7;
ps.mag_heading = random('Uniform', 0,2*pi,n,1);
ps.prob = ones(n,1)*(1/n);

hold on 
plot(lm.x,lm.y,'.','MarkerSize', 10)
axis image
% sdf(gcf,'sj2')
% plot(shp,'FaceAlpha',.1,'EdgeColor','r')
h_ps = scatter(ps.x,ps.y,20,'c','filled','MarkerFaceAlpha',.2);
h_pm = plot(mean(ps.x),mean(ps.y),'ms');
legend('reference point','particle','estimated location','location','best')
hold off
% ylim([5, 25])
% xlim([5, 83])

frame = getframe(gcf);
writeVideo(v,frame);

tM = testdata.Magnetometer(locs,:);

est = zeros(length(tM),3);

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
    if algorithm == 1
        rotatedMag = getHeadingRotatedVector(ps.mag_heading, tM(i,:), rotMat(:,:,i));
        mag_dist = sqrt(sum((rotatedMag-lM(I,:)).^2,2));
    elseif algorithm == 2   % related work: Horizontal & Vertical components, MaLoc, B_h,B_v    
        observed = [vecnorm(lM(I,1:2),2,2), lM(I,3)];
        measured = [vecnorm(tM(i,1:2),2), tM(i,3)];
        mag_dist2 = pdist2(observed, measured,'mahalanobis',nearestSPD(nancov(observed)));
        mag_dist = mag_dist2';
    else
        disp('ILoA.m - wrong method input arg')
        mag_dist = 0;
    end

    if ~all(mag_dist)
        disp('im·pov·er·ish·ment')      %% TODO: what do I have to? not seem starvation, what was it?
        mag_dist = 1/n;
    end
    ps.prob = 1./(mag_dist);
%     in = isinterior(shp,ps.x,ps.y);
% %     in = inShape(shp,[ps.x,ps.y]);
%     ps.prob(~in) = 0;
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
    % --------------------------------------------------------------------
    set(h_ps,'XData',ps.x,'YData',ps.y)                             % ps result
    set(h_pm,'XData',mean(ps.x),'YData',mean(ps.y));
    drawnow     
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);