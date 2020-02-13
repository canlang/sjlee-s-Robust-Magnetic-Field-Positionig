clearvars; close all

intp_candi = [.1,.2,.3,.5,.8,1,1.2];
trj_idx = 1;
device_name = 'S9';
site_name = 'KI-1F';

filename = sprintf('est-result/%s-s%d-%s-errs.mat'...
    ,site_name,trj_idx,device_name);
load(filename)
errs(8:10,:) = [];

MED_mat = cellfun(@mean, errs);
true_pos_thr = 3;
true_pos = MED_mat<true_pos_thr;

precision_rate = sum(true_pos,2)/size(errs,2)*100;
disp('MED (m) /precision (%)')

N = length(intp_candi);
% N = 3;
% set(0,'DefaultAxesColorOrder',brewermap(N,'Dark2'))
set(0,'DefaultAxesColorOrder',brewermap(N,'Blues'))
figure;
hold on

for i=fliplr(1:size(errs,1))
    terr = MED_mat(i,true_pos(i,:));
    tempErrCell = errs(i,true_pos(i,:))';
    all_err = vertcat(tempErrCell{:});
    
    h = cdfplot((all_err));       % draw N=2000
    if intp_candi(i) == .8
        set(h,'LineWidth',2,'Color','r')
    end
    fprintf('%.2f / %d \n', mean(terr), precision_rate(i));
end
hold off

xlim([0 5])
legend(arrayfun(@(x) num2str(x,'\\delta = %1.1f m'), fliplr(intp_candi)...
    ,'UniformOutput', false))

% legend(arrayfun(@(x) num2str(x,'\\delta - %1.1f m'), [.1,.8,1.2]...
%     ,'UniformOutput', false))
legend('Location','best')

set(gcf,'units','points','position',[800,100,800,600])
sdf(gcf,'sj2')
