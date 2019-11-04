clear;close all;clc;
addpath('xyz file operations')
addpath('Multiprod_2009')


target_rawdata_paths = getNameFolds('input_rawdata');
rawdata = load_rawdata(fullfile('input_rawdata',target_rawdata_paths{36}));
% 24,27,30,35, sameline
% 26,28,32,33 left
% 25,29,30,34 right

%% resample for synchronize
rate = 2e-2;
processed_data = resample_rawdata(rawdata,rate);


%% find step point (step) and time labeling
% threshold should be tuned experimentally to match a person's level
minPeakHeight = std(processed_data.acc_norm);

[pks,locs] = findpeaks(processed_data.acc_norm,'MinPeakDistance',...
    .3/rate,'MinPeakHeight',minPeakHeight);   % .3s 이내의 피크는 무시

%% rebuild trajectory
addpath(genpath('madgwick_algorithm_matlab'));
AHRS = MadgwickAHRS('SamplePeriod', rate, 'Beta', 0.1); % sample rate: 2e-2

quaternion = zeros(length(processed_data.Time), 4);
for t = 1:length(processed_data.Time)
    AHRS.Update(processed_data.Gyroscope(t,:), processed_data.Accelerometer(t,:),...
        processed_data.Accelerometer(t,:));	% gyroscope units must be radians
    quaternion(t, :) = AHRS.Quaternion;
end

euler = quatern2euler(quaternConj(quaternion)) * (180/pi);	% use conjugate for sensor frame relative to Earth and convert to degrees.

%% x-y labeling (location)
estloc = zeros(length(locs),2);
% phi,theta,psi: roll pitch yaw (aerospace notation sequence?)
for i=1:length(locs)
% for i=1:length(estloc)
    t_i = locs(i);
    yaw = euler(t_i,3)+90;
    step = 0.7;
    trM = [step*(1*cosd(yaw) - 0*sind(yaw));step*(1*sind(yaw) + 0*cosd(yaw))];
%     trM = [step*(0*cosd(yaw) - 1*sind(yaw));step*(0*sind(yaw) + 1*cosd(yaw))];
    if all(estloc == 0)
        estloc(i,:) = (trM+[0;0])';
    else
        estloc(i,:) = (trM+estloc(i-1,:)')';
    end
end

%% sensor value labeling (magnetic)

A = quatern2rotMat(quaternConj(quaternion(locs,:)));
% A = quatern2rotMat((quaternion(locs,:)));
for i=1:length(A)
%     A2(:,:,i)=inv(A(:,:,i));
    x = 0*pi/180;
    A2(:,:,i)=[cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]*A(:,:,i);
%     disp(A(:,:,i))
end
% ---------- Naive B
B = processed_data.Magnetometer(locs,:);

for i=1:3
    subplot(3,1,i)
    plot_idx = i;
    plot(B(:,plot_idx))
    xlabel('# step');
    ylabel('\muT');
    title(sprintf('%3d degree',0));
    % sdf(gcf,'sj2')
    % set(gcf,'units','points','position',[100,500,1200,600])
    % print -clipboard -dbitmap

    % ---------- Rotation B
    % bsxfun(@times,A,B)
    % C = arrayfun(@(a) a,A);
    C = multiprod(A2,B',[1,2],1);

    hold on
    h = plot(C(plot_idx,:));
    hold off
    % legend('naive','tranform');
end

% set(gcf,'units','points','position',[500,500,1200,600])
% sdf(gcf,'sj2')

% ---------- Inversed Rotation B
% C2 = zeros(size(C));
% for i=1:length(A2)
%     C2(:,i)=A2(:,:,i)\B(i,:)';
% end
% hold on
% h = plot(C2(plot_idx,:));
% hold off
% legend('naive','tranform');

set(gcf,'units','points','position',[1000,500,1000,1000])
sdf(gcf,'sj2')

% video_flag = 0;
% if video_flag
%     v = VideoWriter(strcat('movie_files/','rotation_offset','.mp4'),'MPEG-4');
%     v.FrameRate = 20;
%     v.Quality = 100;
%     open(v);
%     frame = getframe(gcf);
%     writeVideo(v,frame);
% end
% for j=1:360
%     for i=1:length(A)
%     %     A2(:,:,i)=inv(A(:,:,i));
%         x = -j*pi/180;
%         A2(:,:,i)=[cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]/A(:,:,i);    
%     %     disp(A(:,:,i))
%     end
%     C = multiprod(A2,B',[1,2],1);
% %     h.XData = C(1,:);
%     h.YData = C(1,:);
%     title(sprintf('%3d degree',j));
%     drawnow
%     if video_flag
%         frame = getframe(gcf);
%         writeVideo(v,frame);
%     end
% end
% if video_flag close(v);end    

print -clipboard -dbitmap
return
%%
[X,Y,Z]=xyz2grid(estloc(:,1),estloc(:,2),processed_data.Magnetometer(locs,1));
h = imagesc(X(1,:),Y(:,1),Z);
xlabel('x (m)')
ylabel('y (m)')
axis image
axis xy
set(h,'alphadata',~isnan(Z))
% colormap(h,cool(10))
colorbar

%% local functions
% ------------------------------
% local function#1
% ------------------------------
function nameFolds = getNameFolds(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
% nameFolds(~cellfun(@isempty,regexp(nameFolds, '\d{6}'))) = [];
end
