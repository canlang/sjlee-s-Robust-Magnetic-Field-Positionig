clearvars; close all

intp_candi = [.2,.3,.5,.8,1,1.2];
% intp_candi = [.1,.2,.3,.5,.8,1,1.2];

device_name = 'S9';
site_name = 'KI-1F';

trj_idx = 1;
filename = sprintf('est-result/%s-s%d-%s-errs(100).mat'...
    ,site_name,trj_idx,device_name);
load(filename)
errs1 = errs;

trj_idx = 2;
filename = sprintf('est-result/%s-s%d-%s-errs(100).mat'...
    ,site_name,trj_idx,device_name);
load(filename)
errs2 = errs;

trj_idx = 3;
filename = sprintf('est-result/%s-s%d-%s-errs(100).mat'...
    ,site_name,trj_idx,device_name);
load(filename)
errs3 = errs;

errs = [errs1,errs2,errs3];

% errs(8:end,:) = [];

if trj_idx == 1
    MED_mat = cellfun(@mean, errs);
else
    MED_mat = cellfun(@(x) mean(x(3:end)), errs);
end

true_pos_thr = 3;
true_pos = MED_mat<true_pos_thr;

precision_rate = sum(true_pos,2)/size(errs,2)*100;
disp('MED (m) /precision (%)')

N = length(intp_candi);
% N = 3;
% set(0,'DefaultAxesColorOrder',brewermap(N,'Dark2'))
set(0,'DefaultAxesColorOrder',brewermap(N,'Greys'))
figure;
hold on

for i=fliplr(1:size(errs,1))
    terr = MED_mat(i,true_pos(i,:));
%     if trj_idx == 1
%         tempErrCell = errs(i,true_pos(i,:))';
%     else
%         tempErrCell = errs(i,true_pos(i,2:end))';
%     end
    tp_errs = errs(i,true_pos(i,:));
    tp_errs_flt = cellfun(@(x) x(2:end),tp_errs,'UniformOutput',false);
    all_err = vertcat(tp_errs_flt{:});
    
    h = cdfplot((all_err));       % draw N=2000
    if intp_candi(i) == 1
%         set(h,'LineWidth',2,'Color','r')
    end
    fprintf('%.2f / %3.f \n', mean(terr), precision_rate(i));
end
hold off

% xlim([0 10])
legend(arrayfun(@(x) num2str(x,'\\delta = %1.1f m'), fliplr(intp_candi)...
    ,'UniformOutput', false))

% legend(arrayfun(@(x) num2str(x,'\\delta - %1.1f m'), [.1,.8,1.2]...
%     ,'UniformOutput', false))
legend('Location','best')

set(gcf,'units','points','position',[800,100,800,600])
sdf(gcf,'sj2')
