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

%%
for i=["mamam","fafa"]
    disp(i)
end
%%
clearvars;clc;
% a = magic(3)
% 
% hold on
% plot(a(:,1))
% plot(a(:,2))
% plot(a(:,3))
% hold off


x = linspace(0,7);
y = ones(length(x),9);
for i = 1:9
    y(:,i) = sin(x-i/5)';
end
plot(x,y)

line_type = {'--', '-', ':'}; 
set(gca,'LineStyleOrder',line_type)