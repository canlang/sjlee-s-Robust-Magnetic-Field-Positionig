clearvars;close all;
% plotting
intp_candi = [.1,.2,.3,.5,.8,1,1.2];

N = length(intp_candi);
% N = 3;
% set(0,'DefaultAxesColorOrder',brewermap(N,'Dark2'))
set(0,'DefaultAxesColorOrder',brewermap(N,'Blues'))

figure
hold on

for j = fliplr(1:length(intp_candi))
    filename = sprintf('est-result/n1-7f-parfor-%s.mat',num2str(intp_candi(j)));
    load(filename)
    gt_length = 63;
    true_pos_thr = 3;

    % N=2000 case
    conv_case = convIndexes<gt_length/2;
    errMeanMat = cellfun(@mean,errMat);
    true_pos = errMeanMat<true_pos_thr;

    Nc = [1000,2000,3000];
    fprintf('Interpolation granularity - %.1f\n',intp_candi(j));
    for i=1:3
        v_precision = sum(true_pos(i,:))/sum(conv_case(i,:));
        errs = errMat(i,true_pos(i,:));
        err = vertcat(errs{:});
%         if i==3 && ismember(intp_candi(j),intp_candi)
        if i==2     % draw N=2000
%             h = cdfplot(rmoutliers(err));       
            h = cdfplot((err));       % draw N=2000
            if intp_candi(j) == .8
                set(h,'LineWidth',4,'Color','r')
%             else
%                 cdfplot(rmoutliers(err))
            end
        end
        fprintf('(N=%d) %.2f / %2.0f (MED/Precision)\n'...
            ,Nc(i), mean(err), v_precision*100);
    end
    
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