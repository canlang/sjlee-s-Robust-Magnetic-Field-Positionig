function [estloc, c_lambda, c_eta] = getOptTrajectory(euler,step_idx,gt)
estloc = zeros(length(step_idx),2);
c_lambda = zeros(length(step_idx),1);
c_eta = zeros(length(step_idx),1);

yaw = euler(step_idx,3)*pi/180;
rate = 0.8;        % step rate?
[~,turn_idx] = findpeaks(stdfilt(yaw),'MinPeakDistance',2/rate,'MinPeakHeight',0.4);
const_idx = [1; turn_idx; length(step_idx)];           % constraints, checkpoint, waypoints
    
for i=1:length(gt)-1   
    % const_idx(i+1): pair of (next point) the line, 
    % -1: computing before segment item
    if i ~= length(gt)-1        % t_idx: target
        t_idx = const_idx(i):const_idx(i+1)-1;      
    else                        % last line segment have to contain last location
        t_idx = const_idx(i):const_idx(i+1);        
    end
    
    div_step_idx = step_idx(t_idx);
    div_gt = gt(i:i+1,:);

    if i == 1
        start_loc = div_gt(1,:);
    else
        start_loc = estloc(const_idx(i)-1,:);           % 
    end
    x0 = zeros(length(div_step_idx),1);
    
    fun1 = @(x)sum(myfun(euler,div_step_idx,div_gt,x));
%     opt_lambda = fminsearch(fun1,zeros(length(div_step_idx),1));
    opts = optimoptions('fminunc','Algorithm','quasi-newton');
    opt_lambda = fminunc(fun1,x0,opts);
    c_lambda(t_idx,1) = cumsum(opt_lambda);
%     c_lambda(const_idx(i):const_idx(i+1)-1,1) = cumsum(opt_lambda);

    fun2 = @(x)sum(myfun(euler,div_step_idx,div_gt,opt_lambda,x));
%     opt_eta = fminsearch(fun2,x0);
    opt_eta = fminunc(fun2,x0,opts);
    c_eta(t_idx,1) = cumsum(opt_eta);
%     c_eta(const_idx(i):const_idx(i+1)-1,1) = cumsum(opt_eta);

    estloc(t_idx,:) = getTrajectory(euler,div_step_idx,start_loc,opt_lambda,opt_eta);    
end
%%
function a = myfun(e,si,gt,lambda,eta,varargin)
    if nargin == 4 
        loc = getTrajectory(e,si,gt(1,:),lambda);
        [a, ~] = p_poly_dist(loc(:,1),loc(:,2),gt(:,1),gt(:,2));
    elseif nargin == 5
        loc = getTrajectory(e,si,gt(1,:),lambda,eta);
        a = sqrt(sum((loc(end,:)-gt(end,:)).^2,2));
    end
end

end