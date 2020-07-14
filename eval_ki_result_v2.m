clearvars; close all

% intp_candi = [.1];
% intp_candi = [1,1.2];                 % 201
% intp_candi = [.2,.3,.5,.8,1,1.2];
intp_candi = [.1,.2,.3,.5,.8,1,1.2];    % 20, 100, 

device_name = 'S9';
site_name = 'KI-1F';
n_particle = [3000,4000];
method_name = {'ILoA','MaLoc'};
test_trj = 1:3;

n_color = length(intp_candi);
set(0,'DefaultAxesColorOrder',[colormap(cool(n_color));colormap(copper(n_color))]);

re_errMat = zeros(length(n_particle),length(intp_candi),2);     % 2:= compare two algorithm
re_preMat = re_errMat;

hold on
for i=1:length(n_particle)
    for k=1:2
        fileFormat = 'est-result/%s/n%d/%s_traj_%d-%s-errs-rep(100).mat';
        % 201: [1.0, 1.2], 3000,4000

        errCells = {};
        for trj_idx = test_trj
            filename = sprintf(fileFormat,site_name,n_particle(i),method_name{k},trj_idx,device_name);
            load(filename);
            errCells{trj_idx} = errs;   
        end
        %%
        errs = [errCells{test_trj}];       % can choose only specific trajectory (scenario)
        %%
        % errs(8:end,:) = [];       % ??

        if trj_idx == 1
            MED_mat = cellfun(@mean, errs);
        else
            MED_mat = cellfun(@(x) mean(x(3:end)), errs);
        end

        true_positive_err_thr = 9;                   % right?, what is optimal?
        true_pos = MED_mat<true_positive_err_thr;

        precision_rate = sum(true_pos,2)/size(errs,2)*100;
        disp('MED (m) /precision (%)')

        if ~isequal(size(errs,1),length(intp_candi))
            disp('error: wrong setting intpolation candidate.')
            break
        end
        for j=1:size(errs,1)        % intp
            terr = MED_mat(j,true_pos(j,:));
        %     if trj_idx == 1
        %         tempErrCell = errs(i,true_pos(i,:))';
        %     else
        %         tempErrCell = errs(i,true_pos(i,2:end))';
        %     end
            tp_errs = errs(j,true_pos(j,:));
            tp_errs_flt = cellfun(@(x) x(2:end),tp_errs,'UniformOutput',false);
            all_err = vertcat(tp_errs_flt{:});

            if i==1                     % draw N=3000
                h = cdfplot((all_err));       
            end
            if intp_candi(j) == 1
        %         set(h,'LineWidth',2,'Color','r')
            end            
            fprintf('%.2f /%3.f \n', mean(terr), precision_rate(j));
            re_errMat(i,j,k) = mean(terr);
            re_preMat(i,j,k) = precision_rate(j);
        end
    end
end
hold off
% xlim([0 10])
xlabel('Error distance (m)')
iloa_intp_lgd = arrayfun(@(x) num2str(x,'ILoA (\\delta = %1.1f m)'), intp_candi,'UniformOutput', false);
maloc_intp_lgd = arrayfun(@(x) num2str(x,'MaLoc (\\delta = %1.1f m)'), intp_candi,'UniformOutput', false);
lgd = legend([iloa_intp_lgd,maloc_intp_lgd],'Location','best');
lgd.NumColumns = 2;

set(gcf,'units','points','position',[100,100,800,500])
sdf(gcf,'sj6')
%%
figure
n_color = length(n_particle);
set(0,'DefaultAxesColorOrder',[colormap(cool(n_color));colormap(copper(n_color))])

hold on
for j=1:2
%     plot(re_errMat(3,:,i),'s-','MarkerSize',10)
%     plot(re_errMat(2,:,i),'x-.','MarkerSize',10)
%     plot(re_errMat(1,:,i),'d--','MarkerSize',10)
    if j==1
        plot(re_errMat(1,:,j),'s-','MarkerSize',10)
        plot(re_errMat(2,:,j),'x-','MarkerSize',10)
    else
        plot(re_errMat(1,:,j),'s-.','MarkerSize',10)
        plot(re_errMat(2,:,j),'x-.','MarkerSize',10)
    end
end
grid on
ylim([0 3])

ylabel('MED (m)')
xlabel('\delta (m)')
xticks(1:7)
xticklabels(arrayfun(@(x) num2str(x, '%1.1f'),intp_candi,'UniformOutput',false))
iloa_legend = arrayfun(@(x) num2str(x,'ILoA (N = %d)'), n_particle...
    ,'UniformOutput', false);
maloc_legend = arrayfun(@(x) num2str(x,'MaLoc (N = %d)'), n_particle...
    ,'UniformOutput', false);
lgd = legend([iloa_legend,maloc_legend],'Location','best');
lgd.NumColumns = 2;

set(gcf,'units','points','position',[800,100,800,500])
sdf(gcf,'sj6')

%%
figure
n_color = length(intp_candi);

bar_cmap = [cool(n_color);copper(n_color)];
set(0,'DefaultAxesColorOrder',bar_cmap);

% re_preMat2(:,1) = re_preMat(:,5,2);
% re_preMat2(:,2) = re_preMat(:,5,1);
re_preMat2 = [re_preMat(:,:,1),re_preMat(:,:,2)];
bar(re_preMat2,'EdgeColor','none');
ylabel('Precision (%)')
xlabel('Number of particle (x1000)')
xticklabels({'3','4'})
% set(gca,'XTick',3:4)
% title('precision rate at N1')

iloa_legend = arrayfun(@(x) num2str(x,'ILoA (\\delta = %1.1f m)'), intp_candi,'UniformOutput', false);
maloc_legend = arrayfun(@(x) num2str(x,'MaLoc (\\delta = %1.1f m)'), intp_candi,'UniformOutput', false);
lgd = legend([iloa_legend,maloc_legend],'Location','south');
lgd.NumColumns = 2;

set(gcf,'units','points','position',[1600,100,800,500])
sdf(gcf,'sj6')

%%
iloa_err = re_errMat(:,:,1);
mean_iloa = mean(iloa_err(:));
maloc_err = re_errMat(:,:,2);
mean_maloc = mean(maloc_err(:));
fprintf('ILoA:%1.2f, MaLoc:%1.2f, Imp.%2.1f (%%) \n',...
    mean_iloa,mean_maloc,(mean_maloc-mean_iloa)/mean_maloc*100);
% best_maloc-best_iloa
%% from 201 rep.
% re_errMat(1,6:7,1) = [1.39,1.59];
% re_errMat(2,6:7,1) = [1.30,1.62];
% re_errMat(1,6:7,2) = [1.96,2.59];
% re_errMat(2,6:7,2) = [1.87,2.54];