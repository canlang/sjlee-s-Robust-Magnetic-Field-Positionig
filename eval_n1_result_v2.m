% for MaLoc result combine to compare accuracy
clearvars;close all;
% plotting
intp_candi = [.1,.2,.3,.5,.8,1,1.2];
Nc = [1000,2000,3000];

N = length(intp_candi);
% N = 3;
% set(0,'DefaultAxesColorOrder',brewermap(N,'Dark2'))
set(0,'DefaultAxesColorOrder',brewermap(N+3,'*Greys'))

figure
hold on

re_errMat = zeros(length(Nc),length(intp_candi));
re_preMat = re_errMat;

for j = 1:length(intp_candi)
    filename = sprintf('est-result/n1-7f-parfor-%s.mat',num2str(intp_candi(j)));
    load(filename)
    gt_length = 63;
    true_pos_thr = 3;

    % N=2000 case
    conv_case = convIndexes<gt_length/2;
    errMeanMat = cellfun(@mean,errMat);
    true_pos = errMeanMat<true_pos_thr;
 
    fprintf('Interpolation granularity - %.1f\n',intp_candi(j));
    for i=1:length(Nc)
        v_precision = sum(true_pos(i,:))/sum(conv_case(i,:));
        errs = errMat(i,true_pos(i,:));
        err = vertcat(errs{:});
        
        if i==2     % draw N particles, e.g., N=2000 , i==2
%             h = cdfplot(rmoutliers(err));       
            h = cdfplot((err));   
            if intp_candi(j) == .8
                set(h,'LineWidth',4,'Color','r')
            end
        end
%         subplot(1,3,i)
%         hold on
%         h = cdfplot((err));
%         if intp_candi(j) == .8
%             set(h,'LineWidth',4,'Color','r')
%         end
        
        fprintf('(N=%d) %.2f / %2.0f (MED/Precision)\n'...
            ,Nc(i), mean(err), v_precision*100);
        re_errMat(i,j) = mean(err);
        re_preMat(i,j) = mean(v_precision*100);
    end
    
end
hold off

xlim([0 5])
xlabel('Error distance (m)')
legend(arrayfun(@(x) num2str(x,'\\delta = %1.1f m'), intp_candi...
    ,'UniformOutput', false))

% legend(arrayfun(@(x) num2str(x,'\\delta - %1.1f m'), [.1,.8,1.2]...
%     ,'UniformOutput', false))
legend('Location','best')

set(gcf,'units','points','position',[800,100,800,600])
sdf(gcf,'sj3')

%%
close all
% N = length(intp_candi);
N = 4;
set(0,'DefaultAxesColorOrder',brewermap(N,'RdGy'))
% set(0,'DefaultAxesColorOrder',brewermap(N,'Greys'))
figure
hold on
plot(re_errMat(3,:),'s-','MarkerSize',15)
plot(re_errMat(2,:),'x-.','MarkerSize',15)
plot(re_errMat(1,:),'d--','MarkerSize',15)
grid on
ylim([0 .7])

ylabel('MED (m)')
xlabel('\delta (m)')
xticklabels(arrayfun(@(x) num2str(x, '%1.1f'),intp_candi,'UniformOutput',false))
legend('N=3000','N=2000','N=1000','Location','best')

set(gcf,'units','points','position',[800,100,800,600])
sdf(gcf,'sj3')

figure
% hold on
% plot(re_errMat(3,:),'-x','MarkerSize',15)
% plot(re_errMat(2,:),'-+','MarkerSize',15)
% plot(re_errMat(1,:),'-s','MarkerSize',15)
% grid on
% ylim([0 .7])
re_preMat2 = re_preMat';
re_preMat2 = [re_preMat2(:,3),re_preMat2(:,2),re_preMat2(:,1)];
bar(re_preMat2)
set(gcf,'units','points','position',[800,1100,800,600])
sdf(gcf,'sj3')

ylabel('Precision (%)')
xlabel('\delta (m)')
xticklabels(arrayfun(@(x) num2str(x, '%1.1f'),intp_candi,'UniformOutput',false))
legend('N=3000','N=2000','N=1000','Location','best')