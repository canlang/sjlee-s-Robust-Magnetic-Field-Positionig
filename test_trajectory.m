clearvars;close all;
%%
site_name = 'N1-7F';
layout = jsondecode(fileread(sprintf('map/%s.json',site_name)));
opt_lambda = layout.in(:,1);
y = layout.in(:,2);
shp = polyshape(opt_lambda,y);
if ~isempty(layout.out)
    for i = 1:length(layout.out)
        ox = layout.out{i}(:,1);
        oy = layout.out{i}(:,2);
        shp = subtract(shp,polyshape(ox,oy));
    end
end
plot(shp)
xlabel('x (m)');
ylabel('y (m)');
% ylim([-10 50])
axis image
% axis 'auto'
grid on
grid minor
set(gcf,'units','points','position',[300,300,2000,600])

% shapewrite(shp,'main_n17f_roads.shp')

% x1 = [9 22 22]            % replaced manual GUI: drawpolyline
% y1 = [21 21 8.5]
% roi = drawpolyline;  
% gt = roi.Position;
% save('trajectory/gt_traj','gt')
%%
addpath(genpath('p_poly_dist'));
load('trajectory/gt_traj')
load('trajectory/euler_and_step_idx.mat')
% temporal
gt = [gt(1,:);gt(2,:);gt(3,:);gt(2,:);gt(1,:)];
% loc = getOptTrajectory(euler,step_idx,gt);
% hold on
% plot(loc(:,1),loc(:,2),'x-')
%% 
estloc1 = getTrajectory(euler,step_idx,gt(1,:));
[a1, ~] = p_poly_dist(estloc1(:,1),estloc1(:,2),gt(:,1),gt(:,2));
fprintf('rawdata score: %3.1f\n',sum(a1))

estloc2 = zeros(size(estloc1));
estloc3 = zeros(size(estloc1));
%% divde
for i=1:length(gt)-1
    yaw = euler(step_idx,3)*pi/180;
    rate = 0.8;        % step rate?
    [~,turn_idx] = findpeaks(stdfilt(yaw),'MinPeakDistance',2/rate,'MinPeakHeight',0.4);
    const_idx = [1; turn_idx; length(step_idx)];           % constraints, checkpoint, waypoints

    div_step_idx = step_idx(const_idx(i):const_idx(i+1));
    div_gt = gt(i:i+1,:);

    if i == 1
        start_loc = div_gt(1,:);
    else
        start_loc = estloc3(const_idx(i),:);
    end
    
    x0 = zeros(length(div_step_idx),1);
    
    fun1 = @(x)sum(myfun(euler,div_step_idx,div_gt,x));
%     opt_lambda = fminsearch(fun1,x0);
    opts = optimoptions('fminunc','Algorithm','quasi-newton');
    opt_lambda = fminunc(fun1,x0,opts);
    estloc2(const_idx(i):const_idx(i+1),:) = getTrajectory(euler,div_step_idx,start_loc,opt_lambda);   
    
    fun2 = @(x)sum(myfun(euler,div_step_idx,div_gt,opt_lambda,x));
    opt_eta = fminsearch(fun2,x0);
    estloc3(const_idx(i):const_idx(i+1),:) = getTrajectory(euler,div_step_idx,start_loc,opt_lambda,opt_eta);    
  
    loc1 = estloc2(const_idx(i):const_idx(i+1),:);
    loc2 = estloc3(const_idx(i):const_idx(i+1),:);
    
    subplot(2,2,i)
    hold on
    plot(loc1(:,1),loc1(:,2),'x-')
    plot(loc2(:,1),loc2(:,2),'x-')
    axis equal
end
[a2, ~] = p_poly_dist(estloc2(:,1),estloc2(:,2),gt(:,1),gt(:,2));
fprintf('adjusted data score: %3.1f\n',sum(a2))
[a2, ~] = p_poly_dist(estloc3(:,1),estloc3(:,2),gt(:,1),gt(:,2));
fprintf('adjusted data score: %3.1f\n',sum(a2))
%%
hold on
plot(estloc1(:,1),estloc1(:,2),'x-')
% plot(estloc2(:,1),estloc2(:,2),'x-')
plot(estloc3(:,1),estloc3(:,2),'x-')
hold off
axis equal
legend('rawpath','estloc3','location','best');

% sdf(gcf,'sj2')
%%
% close all
% yaw = euler(step_idx,3)*pi/180;
% subplot(311)
% plot(yaw)
% subplot(312)
% plot(unwrap(yaw))
% subplot(313)
% % plot(stdfilt(yaw))
% rate = 0.8;        % step rate?
% findpeaks(stdfilt(yaw),'MinPeakDistance',3/rate,'MinPeakHeight',0.4);
% % [~,turn_idx] = findpeaks(stdfilt(yaw),'MinPeakDistance',1/rate,'MinPeakHeight',0.4);

%%
function a = myfun(e,si,gt,lambda,eta,varargin)
if nargin == 4 
    estloc = getTrajectory(e,si,gt(1,:),lambda);
    [a, ~] = p_poly_dist(estloc(:,1),estloc(:,2),gt(:,1),gt(:,2));
elseif nargin == 5
    estloc = getTrajectory(e,si,gt(1,:),lambda,eta);
    a = sqrt(sum((estloc(end,:)-gt(end,:)).^2,2));
end
end

