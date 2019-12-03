% lb = [5*pi,-20*pi];
% ub = [20*pi,-4*pi];
% 
% opts = optimoptions('ga','PlotFcn',@gaplotbestf);
% 
% 
% rng(1,'twister') % for reproducibility
% IntCon = 1;
% [x,fval,exitflag] = ga(@rastriginsfcn,2,[],[],[],[],...
%     lb,ub,[],IntCon,opts)

rf2 = @(x)rastriginsfcn(x/10);
% rf2 = @(x)rastriginsfcn(x/10); % objective
x0 = [20,30]; % start point away from the minimum
[xf,ff,flf,of] = fminunc(rf2,x0)