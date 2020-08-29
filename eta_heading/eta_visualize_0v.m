clearvars;clc
% iloa = load('mats/iloa_eta_ki.mat','eta');          % ki scenario 3
% maloc = load('mats/maloc_eta_ki.mat','eta');

iloa = load('mats/iloa_success_n1_7f_190228_161738_020-B23A.mat','eta');          % n1
maloc = load('mats/maloc_success_n1_7f_190228_161738_020-B23A.mat','eta');


% https://kr.mathworks.com/matlabcentral/answers/504001-2d-histogram-plot-for-n-x-m-matrix
x = 1:size(iloa.eta,1);
n = 2000;
xRep = repmat(x, 1, n);    


% caz = -69;
% cel = 27;
% caz = -111;
% cel = 50;
caz = -59;
cel = 68;


close all
figure
% subplot(121)
histogram2(xRep(:),iloa.eta(:),[102 100],...
    'Normalization','count','FaceColor','flat','FaceAlpha',0.8);
xlabel('Step (index)')
ylabel('\eta (rad)')
zlabel('Bin Count (\eta)')


set(gca,'YTick',0:pi/2:2*pi) 
set(gca,'YTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'}) 
view(caz,cel)
set(gcf,'units','points','position',[800,500,800,600])
sdf(gcf,'sj2')

% figure
% % subplot(122)
% y = maloc.eta;
% histogram2(xRep(:),maloc.eta(:),[102 100],...
%     'Normalization','count','FaceColor','flat');
% set(gca,'YTick',0:pi/2:2*pi) 
% set(gca,'YTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'}) 
% view(caz,cel)
% % 'FaceColor','flat',
% set(gcf,'units','points','position',[300,200,1800,600])
%%
close all
figure
mean_eta = unwrap(circ_mean(iloa.eta,[],2));
plot(mean_eta,'-')
ylim([0 2*pi])
set(gca,'YTick',0:pi/2:2*pi) 
set(gca,'YTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
% set(gca,'YTick',0:pi/2:pi)
% set(gca,'YTickLabel',{'0','\pi/2','\pi'})
% yline(pi,'r-','y = \pi','LineWidth',3,...
%     'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
yline(pi+0.139,'m--')
% yline(pi+0.139,'m--','Optimal \eta','LineWidth',3,...
%     'LabelHorizontalAlignment','center','LabelVerticalAlignment','bottom');
ylabel('\eta (rad)')
yyaxis right
hold on
plot(abs(mean_eta-(pi+0.139))*180/pi,':')
% plot(abs(mean_eta-(pi))*180/pi,'g-.')
ylim([0 45])
legend('estimated \eta', 'ground truth','heading error','location','best')
ylabel('(degree)')
xlabel('Step (index)')

set(gcf,'units','points','position',[800, 100, 641, 288])
% set(gcf,'units','points','position',[800,100,800,300])
sdf(gcf,'sj')
