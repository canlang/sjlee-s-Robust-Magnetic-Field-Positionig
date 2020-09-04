% Revised from 'm3axis_particle_filter_3d_std_validation_parfor_with_interp
% previous script was used in TMC manuscript
% only N1? test?

clearvars; close all; clc;

% TEST DATASET
data2 = readtable('20171124 MagCoord3axisData.csv');
% data1 = readtable('batch.csv');
% lM = [data1.magnet_x,data1.magnet_y,data1.magnet_z];

attitueCandidate = [0 45 90];
nRepeat = 100;

convIndexes = zeros(length(attitueCandidate),nRepeat,2);    
errMat = cell(length(attitueCandidate),nRepeat,2);
method_name = {'ILoA', 'MaLoc'};

for l=1:1
    algo_idx = l;
    fprintf('%s testing...\n',method_name{algo_idx})
    %%
    for k = 1:length(attitueCandidate)

        attitude_degree=attitueCandidate(k);
        load(sprintf('mats/magmap-n1-7f-step-wpNCE%d.mat',attitude_degree),'map')
        data1.x = map(:,1);data1.y = map(:,2);
        lM = map(:,3:5);
        parfor j = 1:nRepeat
            % ------------------------
            % initialize particle
            n = 2000;

            % 1. @only road
            rand_idx = randi(length(data1.x),n,1);
            ps_x = data1.x(rand_idx);
            ps_y = data1.y(rand_idx);

            % 2. @all area
            % ps_x = random('Uniform', min(data1.x),max(data1.x),n,1);
            % ps_y = random('Uniform', min(data1.y),max(data1.y),n,1);

            % 3. @initial area
            % ps_x = data2.coord_x(1)+random('normal',0,5,n,1);
            % ps_y = data2.coord_y(1)+random('normal',0,5,n,1);
            % ps_y = 15+random('normal',0,5,n,1);

            ps_mag_heading = random('Uniform', 0,2*pi,n,1);
            ps_phy_heading = random('Uniform', 0,2*pi,n,1);
            ps_prob = ones(n,1)*(1/n);
            % ps_stlng = ones(n,1) + random('Uniform', -.1,.1,n,1);
            % ------------------------


            % test data matrix
            tM = [data2.mag_x,data2.mag_y,data2.mag_z];
            % rotate test data 
            x = -pi/2;
            R = [cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1];
            tM = (R*tM')'-25;

            % test result matrix
            est = zeros(length(tM),3);
            err = zeros(length(tM),n);

            for i = 1:length(tM)
                % ================ predict
            %     ps_x = bsxfun(@(x,y) x + cos(y),ps_x,ps_phy_heading);
            %     ps_y = bsxfun(@(x,y) x + sin(y),ps_y,ps_phy_heading);
            %     ps_x = ps_x + cos(ps_mag_heading-pi/2);     % when heading shifted input date
            %     ps_y = ps_y + sin(ps_mag_heading-pi/2);
                ps_x = ps_x + cos(ps_mag_heading);
                ps_y = ps_y + sin(ps_mag_heading);
            %     ps_x = ps_x + ps_stlng.*cos(ps_heading);
            %     ps_y = ps_y + ps_stlng.*sin(ps_heading);

            % ================ update

                % 1. find (geo-locational) nearest learning data
    %             [phy_dist,I] = pdist2([data1.x,data1.y],[ps_x,ps_y],'euclidean','Smallest',1);
                [phy_dist,I] = findNearestLocation([data1.x,data1.y],[ps_x,ps_y]);

                % 2. calculate magnetic distance
    %             R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]),ps_mag_heading,'UniformOutput',false);
    %             rotatedMag = cell2mat(cellfun(@(x)((x*tM(i,:)')'),R,'UniformOutput',false));
                if algo_idx == 1        %TODO: may seem critical bug (not yet applied hotfix)
    %                 rotatedMag = getHeadingRotatedVector(ps_mag_heading, tM(i,:), rotMat(:,:,i));
                    R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]),ps_mag_heading,'UniformOutput',false);
                    rotatedMag = cell2mat(cellfun(@(x)((x*tM(i,:)')'),R,'UniformOutput',false));

                    % TODO: may be more optimizable (DONE?maybe)
                %     mag_dist = diag(pdist2(rotatedMag,lM(I,:),'euclidean'));
                    mag_dist = sqrt(sum((rotatedMag-lM(I,:)).^2,2));
                %     mag_dist = bsxfun(@(x,y) pdist([x;y]), rotatedMag,lM(I,:));
                elseif algo_idx == 2
                    observed = [vecnorm(lM(I,1:2),2,2), lM(I,3)];
                    measured = [vecnorm(tM(i,1:2),2), tM(i,3)];
                    mag_dist2 = pdist2(observed, measured,'mahalanobis',nearestSPD(cov(observed)));
                    mag_dist = mag_dist2';
                else        % Unexpedted (wrong) case
                    mag_dist = 0;
                    disp('wrong feature algoithm input');
                end

                if ~all(mag_dist)
                    disp('Unexpected: all magnetic distance is zero.')
                    mag_dist = 1/n;
                end
                ps_prob = 1./(mag_dist);
                ps_prob(phy_dist>1) = 0;
                if sum(ps_prob) == 0
                    rand_idx = randi(length(data1.x),n,1);
                    ps_x = data1.x(rand_idx);
                    ps_y = data1.y(rand_idx);
                    ps_mag_heading = random('Uniform', 0,2*pi,n,1);
                    ps_phy_heading = random('Uniform', 0,2*pi,n,1);
                    % TODO : have to re-distribute
                    ps_prob = ones(n,1)*(1/n);
                else
                    ps_prob = ps_prob./sum(ps_prob);
                end

                % ================ resample
                resample_idx = randsample(1:n,n,true,ps_prob);
                phy_move_noise_range = 2;
                ps_x = ps_x(resample_idx) + phy_move_noise_range*rand(n,1) - phy_move_noise_range/2;
                ps_y = ps_y(resample_idx) + phy_move_noise_range*rand(n,1) - phy_move_noise_range/2;
                ps_mag_heading = ps_mag_heading(resample_idx)+random('normal',0,.001,n,1);

                est(i,:) = [mean(ps_x),mean(ps_y),mean(abs(angdiff(ps_mag_heading,0)))];        %x,y,heading
                err(i,:) = pdist2([mean(ps_x),mean(ps_y)],[ps_x,ps_y]);
            end

            err_std = std(err,0,2);
            est_err = sqrt(sum((est(:,1:2)-[data2.coord_x (data2.coord_y+.2)]).^2,2));

            cIdx = find(err_std < 2,1);
            if isempty(cIdx)        % not converged case
                cIdx = length(tM);
            end
            convIndexes(k,j,l) = cIdx;
            if cIdx<(length(data2.coord_x)/2)
                errMat{k,j,l} = est_err(cIdx:end);
            end

    %         errMat{k,j} = est_err;


    %         errMat(k,j) = mean(mean(err(63:end,:)));
    %         stdMat(k,j) = std(std(err(63:end,:)));
        end
    end
end
%%
close all;clc
figure
% % axes('ColorOrder',brewermap(3,'RdYlGn'),'NextPlot','replacechildren')
% subplot(1,3,1)
% boxplot(convIndexes',attitueCandidate)
% xlabel('Number of particles N')
% ylabel('Distance to convergence (m)')
% set(gca,'XTickLabelRotation',45)
% 
% subplot(1, 3, [2 3])

% colormap(brewermap([],''))
mean_err = zeros(2,length(attitueCandidate));
n_color = 3;
set(0,'DefaultAxesColorOrder',[colormap(cool(n_color));colormap(copper(n_color))])
hold on
% linS = {'--','-.','-',':','--'};
% linS = {'-','--',':'};
linS = {'-','--'};
for l = 1:1
    for i = 1:length(attitueCandidate)
        h = cdfplot(cell2mat(errMat(i,:,l)'));
        set(h,'lineStyle',linS{l})
    %     set(h,'lineStyle',linS{i},'linewidth',4)
        mean_err(l,i) = mean(cell2mat(errMat(i,:,l)'));
    end
    fprintf('%5s (MED): atti.0 = %.2f ; atti.45 = %.2f ; atti.90 = %.2f\n', method_name{l},...
        mean_err(l,1),mean_err(l,2),mean_err(l,3));
end
% txt1 = strcat(method_name{1},{'(\Deltaatti.='},int2str(attitueCandidate'),'°)');
% txt2 = strcat(method_name{2},{'(\Deltaatti.='},int2str(attitueCandidate'),'°)');

% txt1 = strcat(num2str((1:3)','ILoA Atti.G.#%d '),{' (MED='},num2str(mean_err(1,:)','%.2f)'));
% txt2 = strcat(num2str((1:3)','MaLoc Atti.G.#%d '),{' (MED='},num2str(mean_err(2,:)','%.2f)'));
% txt = {txt1;txt2};

txt1 = strcat(num2str((1:3)','Atti.G.#%d '),{' (MED='},num2str(mean_err(1,:)','%.2f)'));
txt = {txt1};
lgd = legend(vertcat(txt{:}),'location','best');
% lgd.NumColumns = 2;
xlabel('Error distance (m)')

set(gcf,'units','points','position',[200,500,600,400])
sdf(gcf,'sj4')
% tightfig(gcf);


%%
% save_dir = fullfile('exp_mats','N1-7F-CDF');
% if ~exist(save_dir,'dir')
%     mkdir(save_dir)
% end
% 
% if exist('intp_intv','var')
%     mat_filename = sprintf('%s/%s-n1-7f-parfor-%s.mat',save_dir,method_name{algo_idx},num2str(intp_intv));
% else
%     mat_filename = sprintf('%s/%s-n1-7f-parfor-atti_%d.mat',save_dir,method_name{algo_idx},attitude_degree);
% end
% 
% if ~exist(mat_filename,'file')
%     save(mat_filename,'convIndexes','errMat','nParticleCandidate')
% else
%     fprintf('%s file exist, so this is not saved.\n', mat_filename);
% end
return
%% box ploting ILoA & MaLoc: departed to 'boxplotting_ILoA_MaLoc.m'
% clear g
% 
% for i=1:1
%     disp(method_name{i})
%     mat_filename = sprintf('%s/%s-n1-7f-parfor-%s.mat',save_dir,method_name{i},num2str(interp_interval));
%     file_struct = load(mat_filename);
%     fn = fieldnames(file_struct);   
%     % fn{3}: case of the number of particles, fn{2}: error matrix, fn{1}: index of
%     % considring converged    
%     err_mat = file_struct.(fn{2});
%     n_ptcls = file_struct.(fn{3});
%     
%     label_n_particles = {};
%     err = {};
%     label_method = {};
%     for j=1:length(file_struct.(fn{3}))
%         err{j} = cell2mat(err_mat(j,:)')
%         label_n_particles{j} = ones(size(err{j}))*n_ptcls(j)
%         label_method{j} = ones(size(err{j}))*i
% %         label_method{j} = repmat(method_name{i},size(err{j}));
%     end
% 
% end
% % label_n_particles = cell2mat(label_n_particles');
% g(1,1)=gramm('x',label_n_particles,'y',err);
% g(1,1).stat_boxplot();
% g(1,1).set_title('stat_boxplot()');
% figure('Position',[100 100 800 800]);
% g.draw();
%%
% plot(nParticleCandidate,mean(errMat,2),'x-')

boxPlot = @iosr.statistics.boxPlot;
% boxPlot(errMat','style','hierarchy','groupLabels',num2cell(string(nParticleCandidate)))
% boxPlot(errMat','style','hierarchy','groupLabels',nParticleCandidate)
figure('color','w');
h2 = boxPlot(attitueCandidate,errMat', ...
    'boxcolor','auto',...
    'showScatter',true);
%     'showViolin', true, 'boxWidth', 0.025, 'showOutliers', false);
    

    
box on

xlabel('n = # of Particles')
ylabel('Error distance (m)')
sdf(gcf,'sj2')
% set(gca,'color','none')
print -clipboard -dbitmap
% set(gcf,'units','points','position',[800,300,800,800])
% sdf(gcf,'sj2')
% print -clipboard -dbitmap
return
%%
notBoxPlot(errMat(5:end,:)',attitueCandidate(5:end),'style','line','markMedian',true,'interval','SEM')