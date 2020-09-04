% for MaLoc result combine to compare accuracy
clearvars;close all;
% plotting

method_name = {'ILoA', 'MaLoc'};

intp_candi = [.1,.2,.3,.5,.8,1,1.2];
n_particle = [1000,2000,3000];
    
re_errMat = zeros(length(n_particle),length(intp_candi),2);     % 2:= compare two algorithm
re_preMat = re_errMat;

figure
for k = (1:2)
    algo_idx = k;
    n_color = length(intp_candi);
    disp(method_name{k})
%     if k == 1
%         set(groot,'DefaultAxesColorOrder',brewermap(N*2,'*RdBu'))
%         ax = gca;
%         ax.ColorOrderIndex = 1;        
%     else        
%         set(groot,'DefaultAxesColorOrder',brewermap(N*2,'*PuOr'))
%     end   
    set(0,'DefaultAxesColorOrder',[colormap(cool(n_color));colormap(copper(n_color))])
    

    for j = 1:length(intp_candi)
        % the mat file was from the '~~std_validation_parfor_with_interp.m'
        filename = sprintf('exp_mats/N1-7F-CDF/set2/%s-n1-7f-parfor-%s.mat',method_name{algo_idx},num2str(intp_candi(j)));
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

            if i==3   % draw N particles, e.g., N=2000 , i==2
    %             h = cdfplot(rmoutliers(err));                      
                h = cdfplot((err));
                hold on
                if k == 2
                    set(h,'LineStyle','-.');
                end
%                 if intp_candi(j) == .8 && k == 1
%                     set(h,'LineWidth',4,'Color','r')
%                 end
            end
            
%             subplot(1,3,i)
%             hold on
%             h = cdfplot((err));
% %             if intp_candi(j) == .8
% %                 set(h,'LineWidth',4,'Color','r')
% %             end
%             xlabel('Error distance (m)')
%             xlim([0 10])
                        
            fprintf('(N=%d) %.2f / %2.0f\n',n_particle(i), mean(err), v_precision*100); % (MED/Precision)
            re_errMat(i,j,k) = mean(err);
            re_preMat(i,j,k) = mean(v_precision*100);
        end    
    end
end
hold off

xlim([0 10])
xlabel('Error distance (m)')

iloa_legend = arrayfun(@(x) num2str(x,'ILoA (\\delta = %1.1f m)'), intp_candi...
    ,'UniformOutput', false);
maloc_legend = arrayfun(@(x) num2str(x,'MaLoc (\\delta = %1.1f m)'), intp_candi...
    ,'UniformOutput', false);


lgd = legend([iloa_legend,maloc_legend],'Location','best');
lgd.NumColumns = 2;

% legend(arrayfun(@(x) num2str(x,'\\delta - %1.1f m'), [.1,.8,1.2]...
%     ,'UniformOutput', false))

set(gcf,'units','points','position',[800,100,800,500])
% set(gcf,'units','points','position',[800,100,2400,600])
sdf(gcf,'sj6')

%%
% close all
% N = length(intp_candi);

% N = 4;
% set(0,'DefaultAxesColorOrder',brewermap(N,'RdGy'))

n_color = 3;
set(0,'DefaultAxesColorOrder',[colormap(cool(n_color));colormap(copper(n_color))])

figure
hold on
for i=1:2
%     plot(re_errMat(3,:,i),'s-','MarkerSize',10)
%     plot(re_errMat(2,:,i),'x-.','MarkerSize',10)
%     plot(re_errMat(1,:,i),'d--','MarkerSize',10)
    if i==1
        plot(re_errMat(1,:,i),'s-','MarkerSize',10)
        plot(re_errMat(2,:,i),'x-','MarkerSize',10)
        plot(re_errMat(3,:,i),'d-','MarkerSize',10)
    else
        plot(re_errMat(1,:,i),'s-.','MarkerSize',10)
        plot(re_errMat(2,:,i),'x-.','MarkerSize',10)
        plot(re_errMat(3,:,i),'d-.','MarkerSize',10)
    end
end
grid on
ylim([0 2])

ylabel('MED (m)')
xlabel('\delta (m)')
xticklabels(arrayfun(@(x) num2str(x, '%1.1f'),intp_candi,'UniformOutput',false))
iloa_legend = arrayfun(@(x) num2str(x,'ILoA (N = %d)'), n_particle...
    ,'UniformOutput', false);
maloc_legend = arrayfun(@(x) num2str(x,'MaLoc (N = %d)'), n_particle...
    ,'UniformOutput', false);
lgd = legend([iloa_legend,maloc_legend],'Location','best');
lgd.NumColumns = 2;

set(gcf,'units','points','position',[800,100,600,400])
sdf(gcf,'paperwidth')
%%
clear re_preMat2
figure
n_color = 7;

% bar_cmap = flipud([cool(n_color);copper(n_color)]);
% bar_cmap = bar_cmap([1,n_color+1],:);

bar_cmap = [cool(n_color);copper(n_color)];
set(0,'DefaultAxesColorOrder',bar_cmap);
% for i=1:2
%     hold on
%     % plot(re_errMat(3,:),'-x','MarkerSize',15)
%     % plot(re_errMat(2,:),'-+','MarkerSize',15)
%     % plot(re_errMat(1,:),'-s','MarkerSize',15)
%     % grid on
%     % ylim([0 .7])
%     
%     re_preMat2 = re_preMat(:,:,i)';
%     % re_preMat2 = [re_preMat2(:,3),re_preMat2(:,2),re_preMat2(:,1)];     %% why reverse?
%     re_preMat2 = [re_preMat2(:,1),re_preMat2(:,2),re_preMat2(:,3)];    
%     bar(re_preMat2,'FaceAlpha',0.5)
% end
% ylabel('Precision (%)')
% xlabel('\delta (m)')

% bar3([re_preMat(:,:,1)',re_preMat(:,:,2)'])

% re_preMat2(:,1:2:5) = re_preMat(:,:,2)';
% re_preMat2(:,2:2:6) = re_preMat(:,:,1)';
% bar3(re_preMat2,0.5);
% ylabel('Precision (%)')
% ylabel('\delta (m)')
% yticklabels(arrayfun(@(x) num2str(x, '%1.1f'),intp_candi,'UniformOutput',false))
% bar_lgd(1:2:5) = maloc_legend;
% bar_lgd(2:2:6) = iloa_legend;

% re_preMat2(:,1) = re_preMat(:,5,2);
% re_preMat2(:,2) = re_preMat(:,5,1);
re_preMat2 = [re_preMat(:,:,1),re_preMat(:,:,2)];
bar(re_preMat2,'EdgeColor','none');
ylabel('Precision (%)')
xlabel('Number of particle (x1000)')
% title('precision rate at N1')

iloa_legend = arrayfun(@(x) num2str(x,'ILoA (\\delta = %1.1f m)'), intp_candi,'UniformOutput', false);
maloc_legend = arrayfun(@(x) num2str(x,'MaLoc (\\delta = %1.1f m)'), intp_candi,'UniformOutput', false);
lgd = legend([iloa_legend,maloc_legend],'Location','south');
lgd.NumColumns = 2;

set(gcf,'units','points','position',[800,1100,600,400])
sdf(gcf,'paperwidth')


%%
iloa_err = re_errMat(:,:,1);
mean_iloa = mean(iloa_err(:));
maloc_err = re_errMat(:,:,2);
mean_maloc = mean(maloc_err(:));
fprintf('ILoA:%1.2f, MaLoc:%1.2f, Imp.%2.1f (%%) \n',...
    mean_iloa,mean_maloc,(mean_maloc-mean_iloa)/mean_maloc*100);
