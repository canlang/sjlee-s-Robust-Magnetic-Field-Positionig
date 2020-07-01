clear; close all;
data1 = readtable('batch.csv');
data2 = readtable('20171124 MagCoord3axisData.csv');
lM = [data1.magnet_x,data1.magnet_y,data1.magnet_z];

% Video save
% video_flag = chooseVideoSaveOrNot;
video_flag = false;


A = imread('map/N1-7F.png','BackgroundColor',[1 1 1]);

xWorldLimits = [-1 1650/20];
yWorldLimits = [-1 660/20];
RA = imref2d(size(A),xWorldLimits,yWorldLimits);
% imshow(flipud(A),RA);
axis image tight;


% draw learning data
hold on
% plot(data1.x,data1.y,'.')
scatter(data1.x,data1.y,'o','MarkerEdgeColor','b',...
    'MarkerFaceColor',[.49 1 .63],'MarkerFaceAlpha',.6)
xlabel('(m)')
ylabel('y (m)')
% for save eps

set(gcf,'units','points','position',[700,500,1000,350])
% tightfig
legend('reference point')
% legend('reference point','location','best')   
sdf(gcf,'sj4')
print -depsc2 eps/env_setting.eps

% axis equal
% xlim([8 83])
% ylim([0 30])
% set(gcf,'units','points','position',[700,500,1500,700])

%%
% initialize particle
n = 2000;
% 1. only road
rand_idx = randi(length(data1.x),n,1);
ps.x = data1.x(rand_idx);
ps.y = data1.y(rand_idx);

% 2. all area
% ps.x = random('Uniform', min(data1.x),max(data1.x),n,1);
% ps.y = random('Uniform', min(data1.y),max(data1.y),n,1);

% 3. initial area
% ps.x = data2.coord_x(1)+random('normal',0,5,n,1);
% ps.y = data2.coord_y(1)+random('normal',0,5,n,1);
% ps.y = 15+random('normal',0,5,n,1);

ps.mag_heading = random('Uniform', 0,2*pi,n,1);
ps.phy_heading = random('Uniform', 0,2*pi,n,1);
ps.prob = ones(n,1)*(1/n);
% ps.stlng = ones(n,1) + random('Uniform', -.1,.1,n,1);


% draw particle
hold on 
% h_ps = plot(ps.x,ps.y,'.','MarkerSize',8);
h_ps = scatter(ps.x,ps.y,20,'c','filled','MarkerFaceAlpha',.2);
h_pm = plot(mean(ps.x),mean(ps.y),'ms');
h_gt = plot(data2.coord_x(1),data2.coord_y(1),'s','MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
hold off


% test data matrix
tM = [data2.mag_x,data2.mag_y,data2.mag_z];
% rotate test data 
x = -pi/2;

R = [cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1];
% tM = (R*tM')';
tM = (R*tM')'-25;
% 25 offset meaning TestData shift: because magnetometer's calibration not matched

% test result matrix
est = zeros(length(tM),3);
err = zeros(length(tM),n);


if video_flag
    v = VideoWriter('m3axis_2d_pf_nonshiftedinput.mp4','MPEG-4');
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
    ps.x = ps.x + cos(ps.mag_heading);
    ps.y = ps.y + sin(ps.mag_heading);
%     ps.x = ps.x + ps.stlng.*cos(ps.heading);
%     ps.y = ps.y + ps.stlng.*sin(ps.heading);
    
    % ================ UPDATE
    
    % 1. find (geo-locational) nearest learning data
    [phy_dist,I] = findNearestLocation_mex([data1.x,data1.y],[ps.x,ps.y]);
       
    % 2. calculate magnetic distance
    % Converted 'getHeadingRotatedVector' function
%     R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]),ps.mag_heading,'UniformOutput',false);
%     rotatedMag = cell2mat(cellfun(@(x)((x*tM(i,:)')'),R,'UniformOutput',false));
    rotatedMag = getHeadingRotatedVector(ps.mag_heading, tM(i,:));
%     if ~isequal(rotatedMag,rotatedMag2)
%         disp('error convert')
%     end

    % TODO: may be more optimizable (DONE?maybe)
%     mag_dist = diag(pdist2(rotatedMag,lM(I,:),'euclidean'));
    mag_dist = sqrt(sum((rotatedMag-lM(I,:)).^2,2));
%     mag_dist = bsxfun(@(x,y) pdist([x;y]), rotatedMag,lM(I,:));

    % related work: Horizontal & Vertical components, MaLoc, B_h,B_v
%     mag_dist2 = mahal([vecnorm(tM(i,1:2),2), tM(i,3)], [vecnorm(lM(I,1:2),2,2), lM(I,3)]);
%     bsxfun(@(x) mahal([vecnorm(tM(i,1:2),2), tM(i,3)],x), );
    mag_dist2 = pdist2([vecnorm(lM(I,1:2),2,2), lM(I,3)], [vecnorm(tM(i,1:2),2), tM(i,3)],'mahalanobis');
    mag_dist = mag_dist2';

    if ~all(mag_dist)
        break
    end
    ps.prob = 1./(mag_dist);
    ps.prob(phy_dist>1) = 0;
    if sum(ps.prob) == 0
        rand_idx = randi(length(data1.x),n,1);
        ps.x = data1.x(rand_idx);
        ps.y = data1.y(rand_idx);
        ps.mag_heading = random('Uniform', 0,2*pi,n,1);
        ps.phy_heading = random('Uniform', 0,2*pi,n,1);
        ps.prob = ones(n,1)*(1/n);
    else
        ps.prob = ps.prob./sum(ps.prob);
    end
    
    
    % ================ resample
    resample_idx = randsample(1:n,n,true,ps.prob);
    phy_move_noise_range = 2;
%     ps.x = ps.x(resample_idx) + phy_move_noise_range*rand(n,1) - phy_move_noise_range/2;
%     ps.y = ps.y(resample_idx) + phy_move_noise_range*rand(n,1) - phy_move_noise_range/2;
    ps.x = ps.x(resample_idx) + random('normal',0,.5,n,1);
    ps.y = ps.y(resample_idx) + random('normal',0,.5,n,1);
    ps.mag_heading = ps.mag_heading(resample_idx)+random('normal',0,.001,n,1);
%     ps.phy_heading = ps.phy_heading(resample_idx)+random('normal',0,.001,n,1);

%     ps.heading = ps.heading(resample_idx)+random('Uniform', -pi/10,pi/10,n,1);
%     ps.stlng = ps.stlng(resample_idx) + random('normal',0,1,n,1);
%     ps.prob = init_prob;
   
    set(h_ps,'XData',ps.x,'YData',ps.y)                             % ps result
    set(h_pm,'XData',mean(ps.x),'YData',mean(ps.y));
    set(h_gt,'XData',data2.coord_x(i),'YData',data2.coord_y(i))     % ground truth
    drawnow limitrate
    
%     est(i,:) = [mean(ps.x),mean(ps.y),mean(angdiff(ps.mag_heading,0))];
    est(i,:) = [mean(ps.x),mean(ps.y),mean(abs(angdiff(ps.mag_heading,0)))];        %x,y,heading
    err(i,:) = pdist2([data2.coord_x(i),data2.coord_y(i)],[ps.x,ps.y]);
%     break
%     pause(.1)
    
    if video_flag
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
    ps_hist(i) = ps;
end
if video_flag; close(v); end
disp ('current: MaLoc result')
%%
figure
[lineh, bandsh] = fanChart(1:size(err,1),err, 'mean', 10:10:90, ...
    'alpha', .2, 'colormap', {'shadesOfColor', [0 0 .8]});
txt = strcat({'Pct'}, cellstr(int2str((20:20:80)')));
legend([bandsh;lineh], [txt;{'Mean'}])
% clickableLegend([lineh;bandsh], [{'Mean'};txt])

return 
%% visualize magnetic heading (temporal, not completed)
% mean_mag_heading = [];
% for i = 1:length(ps_hist)
%     ps = ps_hist(i);
%     mean_mag_heading(end+1) = mean(ps.mag_heading);
% end
% PC1 = 1.*cos(mean_mag_heading)+pi/2);
% PC2 = 1.*sin(mean_mag_heading)+pi/2);
% s_plot = scatter3(ps.x,PC1(:),PC2(:),100,ps.mag_heading,'filled');

%% error ellipse using "gramm"
g(1,1)=gramm('x',ps.x,'y',ps.y);
g(1,1).geom_point();

g(1,1).set_point_options('base_size',2);
g(1,1).stat_ellipse('type','95percentile','geom','area','patch_opts',{'FaceAlpha',0.1,'LineWidth',2});
g(1,1).set_title('stat_ellispe()');

figure('Position',[100 100 800 600])
g.draw();


%%
figure
A = imread('N1-7F.png','BackgroundColor',[1 1 1]);

xWorldLimits = [-1 1650/20];
yWorldLimits = [-1 660/20];
RA = imref2d(size(A),xWorldLimits,yWorldLimits);
imshow(flipud(A),RA);
axis xy;

hold on
plot(est(:,1), est(:,2),'*r')

%%
figure
subplot(211)
histogram(wrapTo2Pi(ps.phy_heading))
xlabel('Physical Heading')
xlim([0 2*pi])
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'}) 
subplot(212)
% stem(wrapTo2Pi(ps.mag_heading),ones(n,1))
histogram(wrapTo2Pi(ps.mag_heading),n/2)
xlabel('Magnetic Heading')
xlim([0 2*pi])
set(gca,'XTick',0:pi/2:2*pi) 
set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'})

set(gcf,'units','points','position',[800,300,800,500])
sdf(gcf,'sj')


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