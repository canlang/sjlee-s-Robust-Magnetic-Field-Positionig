% for MaLoc result combine to compare accuracy
clearvars;close all;clc;

atti_candi = [0,45,85];
n_particle = [1000,2000,3000];
    
re_errMat = zeros(length(n_particle),length(atti_candi));     % 2:= compare two algorithm
re_preMat = re_errMat;

figure    
% n_color = length(atti_candi);
% set(0,'DefaultAxesColorOrder',[colormap(cool(n_color));colormap(copper(n_color))])

for j = 1:length(atti_candi)
    filename = sprintf('exp_mats/N1-7F-CDF/atti-50rep-3/MaLoc-n1-7f-parfor-atti_%d.mat',atti_candi(j));
    load(filename)
    gt_length = 63;
    true_pos_thr = 3;

    % N=2000 case
    conv_case = convIndexes<gt_length/2;
    errMeanMat = cellfun(@mean,errMat);
    true_pos = errMeanMat<true_pos_thr;

%         fprintf('Interpolation granularity - %.1f\n',intp_candi(j));
    for i=1:length(n_particle)
        v_precision = sum(true_pos(i,:))/sum(conv_case(i,:));
        errs = errMat(i,true_pos(i,:));
        err = vertcat(errs{:});

        if i==2   % draw N particles, e.g., N=2000 , i==2
            h = cdfplot(rmoutliers(err));                      
%             h = cdfplot((err));
            hold on
        end
        fprintf('(N=%d) %.2f / %2.0f\n',n_particle(i), mean(err), v_precision*100); % (MED/Precision)
        re_errMat(i,j) = mean(err);
        re_preMat(i,j) = mean(v_precision*100);
    end    
end

hold off

xlim([0 10])
xlabel('Error distance (m)')

iloa_legend = arrayfun(@(x) num2str(x,'Attitude = %d°'), atti_candi...
    ,'UniformOutput', false);

lgd = legend(iloa_legend,'Location','best');
% lgd.NumColumns = 2;

set(gcf,'units','points','position',[800,100,800,500])
% set(gcf,'units','points','position',[800,100,2400,600])
sdf(gcf,'sj6')