
t_idx = 1;

iloa_filename = sprintf('iloa_success_ki_tr%d',t_idx);
maloc_filename = sprintf('maloc_success_ki_tr%d',t_idx);
iloa = load(sprintf('mats/%s.mat',iloa_filename),'eta');
maloc = load(sprintf('mats/%s.mat',maloc_filename),'eta');

iloa.eta = iloa.eta(1:70,:);
maloc.eta = maloc.eta(1:70,:);

x = 1:size(iloa.eta,1);
n = 2000;
xRep = repmat(x, 1, n);    

%%
close all
% subplot(211)
histogram2(xRep(:),iloa.eta(:),[size(iloa.eta,1), 100],...
    'DisplayStyle','tile','Normalization','probability','ShowEmptyBins','on');
% ,'ShowEmptyBins','on'
ylim([0, 2*pi])
set(gca,'YTick',0:pi/2:2*pi) 
set(gca,'YTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'}) 
ylabel('\eta (rad)')
xlabel('Step (index)')
c = colorbar;
% ylabel( c, 'Probability' );
% caxis([0, 1e-3]);
set(gca,'ColorScale','log')
caxis([1e-5, 1e-2]);

set(gcf,'units','points','position',[300,600,700,500])
% set(gcf,'units','points','position',[300,600,900,300])
sdf(gcf,'paper_f150')
% tightfig
% print(gcf, '-depsc2', sprintf('eps/%s.eps',iloa_filename))
print(gcf, '-djpeg', sprintf('../eps/%s.jpg',iloa_filename),'-r300')

figure
% subplot(212)
histogram2(xRep(:),maloc.eta(:),[size(maloc.eta,1), 100],...
    'DisplayStyle','tile','Normalization','probability','ShowEmptyBins','on');
ylim([0, 2*pi])
set(gca,'YTick',0:pi/2:2*pi) 
set(gca,'YTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'}) 

xlabel('Step (index)')
ylabel('\eta (rad)')
c = colorbar;
% ylabel( c, 'Probability' );
set(gca,'ColorScale','log')
caxis([1e-5, 1e-2]);

set(gcf,'units','points','position',[1100,600,700,500])
% sdf(gcf,'sj5')
sdf(gcf,'paper_f150')
% tightfig
% print(gcf, '-depsc2', sprintf('../eps/%s.eps',maloc_filename))
print(gcf, '-djpeg', sprintf('../eps/%s.jpg',maloc_filename),'-r300')


