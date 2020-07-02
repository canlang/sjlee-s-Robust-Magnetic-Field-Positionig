%% box ploting ILoA & MaLoc
clearvars;

method_name = {'ILoA', 'MaLoc'};
save_dir = fullfile('exp_mats','N1-7F-CDF');
interp_interval = .8;

for i=1:2
    disp(method_name{i})
    mat_filename = sprintf('%s/%s-n1-7f-parfor-%s.mat',save_dir,method_name{i},num2str(interp_interval));
    file_struct = load(mat_filename);
    fn = fieldnames(file_struct);   
    % fn{3}: case of the number of particles, fn{2}: error matrix, fn{1}: index of
    % considring converged    
    convidx{i} = file_struct.(fn{1});
    err_mat{i} = file_struct.(fn{2});
    n_ptcls{i} = file_struct.(fn{3});
    
    % CASE #0??
%     label_n_particles = {};
%     err = {};
%     label_method = {};    
%     for j=1:length(file_struct.(fn{3}))
%         err{j} = cell2mat(err_mat(j,:)');
% %         label_n_particles{j} = ones(size(err{j}))*n_ptcls(j)
%         label_n_particles{j} = repmat({num2str(n_ptcls(j))},size(err{j}));
% %         label_method{j} = ones(size(err{j}))*i
%         label_method{j} = repmat({method_name{i}},size(err{j}));
%     end
    label_n_ptcl{i} = repmat({'1000';'2000';'3000'},1,length(convidx{i}));
    label_method{i} = repmat(method_name(i), size(convidx{i}));
end
%% case#0
% label_n_particles = cat(1, label_n_particles{:});
% label_method = cat(1, label_method{:});
% err = cat(1,err{:});
% 
% % label_n_particles = cell2mat(label_n_particles');
% g(1,1)=gramm('x',label_n_particles,'y',err,'color',label_method);
%% case#converge index
convidx = {cat(1, convidx{:})};
label_n_ptcl = {cat(1, label_n_ptcl{:})};
label_method = {cat(1, label_method{:})};

clear g;close all
g(1,1)=gramm('x',label_n_ptcl{:}(:),'y',convidx{:}(:),'color',label_method{:}(:));
% g(1,1).stat_boxplot();
% g(1,1).set_title('stat_boxplot()');
% figure('Position',[100 100 800 400]);
% g.draw();
%%
g(1,2)=copy(g(1));
g(1,3)=copy(g(1));
g(2,1)=copy(g(1));
g(2,2)=copy(g(1));

%Raw data as scatter plot
g(1,1).geom_point();
g(1,1).set_title('geom_point()');

%Jittered scatter plot
g(1,2).geom_jitter('width',0.4,'height',0);
g(1,2).set_title('geom_jitter()');

%Averages with confidence interval
g(1,3).stat_summary('geom',{'bar','black_errorbar'});
g(1,3).set_title('stat_summary()');

%Boxplots
g(2,1).stat_boxplot();
g(2,1).set_title('stat_boxplot()');

%Violin plots
g(2,2).stat_violin('fill','transparent');
g(2,2).set_title('stat_violin()');
% figure('Position',[100 100 800 800]);
g.draw();