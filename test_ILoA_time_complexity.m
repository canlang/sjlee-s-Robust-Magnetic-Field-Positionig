clearvars;close all

maps = {};
site_name = 'KI-1F';
mypath = 'mats';        % magnetic map directory
intp = [.1, .2, .3, .5, .8, 1.0, 1.2];
map_size = zeros(1, length(intp));
for i=1:length(intp)
    filename = sprintf('magmap-%s-%.1fa.mat',site_name,intp(i));
    if exist(fullfile(mypath,filename), 'file') == 2
        load(fullfile(mypath,filename), 'map')    
        maps{i} = map;
        map_size(i) = length(map);
    end  
end

n_color = length(intp);
% set(0,'DefaultAxesColorOrder',colormap(pink(n_color)));
% set(0,'DefaultAxesColorOrder',colormap(brewermap(n_color,'PuBu')));
% set(0,'DefaultAxesColorOrder',colormap(brewermap(n_color,'YlOrRd')));
% colormap(brewermap(n_color,'YlOrRd'))

% subplot(121)
b = semilogy(intp,map_size,'o','markerfacecolor','b');
grid on
fit = fit(intp',map_size','power1');
hold on
plot(fit,intp,map_size)
strValues = strtrim(cellstr(num2str(map_size(:),'%d')));
text(intp,map_size,strValues,'VerticalAlignment','bottom');
hold off
uistack(b,'top')

xticks(intp)
%%% setting for bar plotting 
% xticklabels(intp)
% color change
% bar_cmap = brewermap(n_color,'YlOrRd');
% b.FaceColor = 'flat';
% for i=1:length(intp)    
%     b.CData(i,:) = bar_cmap(i,:);
% end
%%%

xlabel('\delta (m)')
ylabel('Number of references')

% subplot(122)
% p = pie(map_size);
% labels = arrayfun(@(x) num2str(x,'#M (\\delta = %1.1f)'),intp,'UniformOutput', false);
% lgd = legend(labels,'position',[0.5525 0.6779 0.1446 0.2450]);
% xlabel('The search space size')
% 
% pText = findobj(p,'Type','text');
% percentValues = get(pText,'String'); 
% 
% for i=5:length(intp)
%     pText(i).String = '';
% end

set(gcf,'units','points','position',[300,100,800,600])
sdf(gcf,'sj2')
% tightfig;