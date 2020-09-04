function estloc = getTrajectory(euler,step_idx,start_loc,lambda,eta,varargin)
% eular phi,theta,psi: roll pitch yaw (aerospace sequence?)

step = 0.7*ones(length(step_idx),1);

estloc = zeros(length(step_idx)+1,2);
if nargin >= 3
    estloc(1,:) = start_loc;
end
if nargin >= 4
    % 1. first version
%     euler(step_idx,3) = euler(step_idx,3) + lambda;
    % 2. weight of yaw chaining, because i assumed: weight event is not independent
    for i=1:length(lambda)        
        euler(step_idx(i:end),3) = euler(step_idx(i:end),3) + lambda(i);
    end
end
if nargin == 5
    for i=1:length(eta)        
        step(i:end) = step(i:end) + eta(i);
    end
end

% TODO: seemed it possible without this loop.
for i=2:length(step_idx)+1
    t_i = step_idx(i-1);
    yaw = euler(t_i,3);
    
    % Transition matrix
    trM = [step(i-1)*(1*cosd(yaw) - 0*sind(yaw));step(i-1)*(1*sind(yaw) + 0*cosd(yaw))];
%     trM = [step*(0*cosd(yaw) - 1*sind(yaw));step*(0*sind(yaw) + 1*cosd(yaw))];
    estloc(i,:) = (trM+estloc(i-1,:)')';
%     if i==1
%         estloc(i,:) = trM';
%     else
%         estloc(i,:) = (trM+estloc(i-1,:)')';
%     end
end
estloc(1,:) = [];