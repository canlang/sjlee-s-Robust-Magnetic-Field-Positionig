clearvars;close all;

%% view the shape of trajectories
% close all;
N = 4;
set(0,'DefaultAxesColorOrder',brewermap(N,'Set1'))
% for i = 1:3
%     test_seq_i = i;
%     a = jsondecode(fileread(sprintf('est-result/ki-gt-s%d.json',test_seq_i)));
% %     subplot(3,1,i)
%     plot(a(:,1),a(:,2),'-')
%     hold on
% %     axis square
%     daspect([1 1 1])
% end
% sdf(gcf,'sj2')
% set(gca,'visible','off')
% set(gcf,'units','points','position',[784   723   776   615])
%%
clearvars;close all

for i = 1:3
    test_seq_i = i;
    gt_filename = sprintf('est-result/ki-gt-s%d.json',test_seq_i);
    gt = jsondecode(fileread(gt_filename));
%     figure
    plot(gt(:,1),gt(:,2),'-')
    hold on

    %% ------------ iphone
    est_mat_name = sprintf('est-result/ki-1f-s%d-iphone.mat',test_seq_i);
    d = load(est_mat_name);
    d.converge_idx = d.converge_idx+4;
    converged_est = [d.est(d.converge_idx:end,1), d.est(d.converge_idx:end,2)];
    plot(converged_est(:,1),converged_est(:,2),':')
    % plot(converged_est(:,1),converged_est(:,2),'*')

%     err_dist = arrayfun(@(x,y) p_poly_dist(x,y,gt(:,1)',gt(:,2)')...
%         ,converged_est(:,1),converged_est(:,2));
%     err = abs(err_dist);
%     % ecdf(err)
%     disp(mean(err))

    %% ------------ s9
    est_mat_name = sprintf('est-result/ki-1f-s%d-s9.mat',test_seq_i);
    d = load(est_mat_name);
    converged_est = [d.est(d.converge_idx:end,1), d.est(d.converge_idx:end,2)];
    plot(converged_est(:,1),converged_est(:,2),'--')
    % plot(converged_est(:,1),converged_est(:,2),'o')

%     err_dist = arrayfun(@(x,y) p_poly_dist(x,y,gt(:,1)',gt(:,2)')...
%         , converged_est(:,1),converged_est(:,2));
%     err = abs(err_dist);
%     % ecdf(err)
%     disp(mean(err))

    %% ------------ mate20
    est_mat_name = sprintf('est-result/ki-1f-s%d-mate20.mat',test_seq_i);
    d = load(est_mat_name);
    d.converge_idx = d.converge_idx+10;
    converged_est = [d.est(d.converge_idx:end,1), d.est(d.converge_idx:end,2)];
    plot(converged_est(:,1),converged_est(:,2),'-.')
    % plot(converged_est(:,1),converged_est(:,2),'+')

%     err_dist = arrayfun(@(x,y) p_poly_dist(x,y,gt(:,1)',gt(:,2)')...
%         , converged_est(:,1),converged_est(:,2));
%     err = abs(err_dist);
%     % ecdf(err)
%     disp(mean(err))
%     hold off
    

% d.turn_locs = d.turn_locs;
% plot(d.est(d.turn_locs,1),d.est(d.turn_locs,2),'s')
% 
% tpoints = [d.est(d.turn_locs,1),d.est(d.turn_locs,2);d.est(end,1),d.est(end,2)];
% diag(pdist2(tpoints,gt(2:end,:)))
end
%%
axis equal
% legend('ground truth','iphone','s9','mate20')
legend('Ground truth','Test#1','Test#2','Test#3','location','best')
set(gcf,'units','points','position',[800,600,700,700])
sdf(gcf,'sj2')