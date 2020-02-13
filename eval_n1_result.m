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
        if i==3
            if intp_candi(j) == .8
                h = cdfplot(rmoutliers(err));       % draw N=2000
                set(h,'LineWidth',2,'Color','r')
            else
                cdfplot(rmoutliers(err))       % draw N=2000
            end
        end
        fprintf('(N=%d) %.2f / %2.1f (MED/Precision)\n'...
            ,Nc(i), mean(err), v_precision*100);
    end
    
end
hold off

legend(arrayfun(@(x) num2str(x,'\\delta - %1.1f m'), fliplr(intp_candi)...
    ,'UniformOutput', false))

% legend(arrayfun(@(x) num2str(x,'\\delta - %1.1f m'), [.1,.8,1.2]...
%     ,'UniformOutput', false))
% legend('Location','best')


sdf(gcf,'sj2')