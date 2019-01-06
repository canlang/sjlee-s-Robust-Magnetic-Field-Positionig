clear; close all;
data1 = readtable('batch.csv');
data2 = readtable('20171124 MagCoord3axisData.csv');
lM = [data1.magnet_x,data1.magnet_y,data1.magnet_z];
lM = vecnorm(lM,2,2);
% draw learning data
plot(data1.x,data1.y,'x')
axis equal
xlim([8 83])
ylim([0 30])
set(gcf,'units','points','position',[800,500,800,500])
%%
% initialize particle
n = 100;
rand_idx = randi(length(data1.x),n,1);
ps.x = data1.x(rand_idx);
ps.y = data1.y(rand_idx);

% ps.x = random('Uniform', min(data1.x),max(data1.x),n,1);
% ps.y = random('Uniform', min(data1.y),max(data1.y),n,1);

% ps.x = data2.coord_x(1)+random('normal',0,5,n,1);
% ps.y = data2.coord_y(1)+random('normal',0,5,n,1);
% ps.y = 15+random('normal',0,5,n,1);

ps.phy_heading = random('Uniform', 0,2*pi,n,1);
% ps.phy_heading = random('Uniform', -2/pi,2/pi,n,1);
ps.prob = ones(n,1)*(1/n);
% ps.stlng = ones(n,1) + random('Uniform', -.1,.1,n,1);

% draw particle
hold on 
h_ps = plot(ps.x,ps.y,'.');
h_gt = plot(data2.coord_x(1),data2.coord_y(1),'s','MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
hold off

% test data matrix
tM = [data2.mag_x,data2.mag_y,data2.mag_z];
% rotate test data 
x = -pi/2;
R = [cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1];
tM = (R*tM')'-25;
% make l2Norm
tM = vecnorm(tM,2,2);

% Video save
video_flag = false;
if video_flag
    v = VideoWriter('m1axis_2d_pf_shifted.mp4','MPEG-4');
    v.FrameRate = 10;
    v.Quality = 100;
    open(v);
    frame = getframe(gcf);
    writeVideo(v,frame);
end

for i = 1:length(tM)
    % ================ predict
%     ps.x = bsxfun(@(x,y) x + cos(y),ps.x,ps.heading);
%     ps.y = bsxfun(@(x,y) x + sin(y),ps.y,ps.heading);
    ps.x = ps.x + cos(ps.phy_heading);
    ps.y = ps.y + sin(ps.phy_heading);
%     ps.x = ps.x + ps.stlng.*cos(ps.heading);
%     ps.y = ps.y + ps.stlng.*sin(ps.heading);
    % ================ update
    % find nearest learning data
    [phy_dist,I] = pdist2([data1.x,data1.y],[ps.x,ps.y],'euclidean','Smallest',1);
    % calculate magnetic distance    
    mag_dist = abs(tM(i)-lM(I));

    if all(mag_dist) == 0
        break
    end
    ps.prob = 1./mag_dist;
    
    ZPC = phy_dist>1;%zero_particle_condition
    if ~all(ZPC)
        ps.prob(phy_dist>1) = 0;
        ps.prob = ps.prob./sum(ps.prob);
    else
        ps.prob = ones(n,1)*(1/n);
    end
    
    % ================ resample
    resample_idx = randsample(1:n,n,true,ps.prob);
    phy_move_noise_range = 1;
    ps.x = ps.x(resample_idx) + phy_move_noise_range*rand(n,1) - phy_move_noise_range/2;
    ps.y = ps.y(resample_idx) + phy_move_noise_range*rand(n,1) - phy_move_noise_range/2;
    ps.phy_heading = ps.phy_heading(resample_idx)+random('normal',0,.001,n,1);
%     ps.stlng = ps.stlng(resample_idx) + random('normal',0,1,n,1);
   
    set(h_ps,'XData',ps.x,'YData',ps.y)
    set(h_gt,'XData',data2.coord_x(i),'YData',data2.coord_y(i))
    drawnow
%     pause(.1)
    if video_flag
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end
if video_flag close(v);end

%%
% figure
% subplot(211)
% hist(ps.phy_heading)
% set(gca,'XTick',0:pi/2:2*pi) 
% set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'}) 
% subplot(212)
% hist(ps.mag_heading)
% set(gca,'XTick',0:pi/2:2*pi) 
% set(gca,'XTickLabel',{'0','pi/2','pi','3pi/2','2pi'}) 