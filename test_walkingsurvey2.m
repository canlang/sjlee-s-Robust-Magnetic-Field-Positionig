clear all

%% N1 7F 
% path_coord{1}.x = [9,70.7];
% path_coord{1}.y = 21;
% 
% path_coord{2}.x = [12,80];
% path_coord{2}.y = 8;
% 
% path_coord{3}.x = 22;
% path_coord{3}.y = [21,8];
% 
% path_coord{4}.x = 38.3;
% path_coord{4}.y = [21,8];
% 
% path_coord{5}.x = 54;
% path_coord{5}.y = [21,8];
% 
% path_coord{6}.x = 71;
% path_coord{6}.y = [21,8];
% 
% path_coord{7}.x = [9,70.7];
% path_coord{7}.y = 21.6;
% 
% path_coord{8}.x = [9,70.7];
% path_coord{8}.y = 20.4;
% 
% path_coord{9}.x = [12,80];
% path_coord{9}.y = 8.6;
% 
% path_coord{10}.x = [12,80];
% path_coord{10}.y = 7.4;

% t_input_idx = 43;       % start index for radio-map rawdata.
% radiomap_index = [44:(44+5),57:(57+3)];     % east-direction
% radiomap_index = [64];          % south-direction

%% N1 2F (grid open space)
% path_coord = cell(1,9);
% for i = 1:9
%     path_coord{i}.x = [38,56];
%     path_coord{i}.y = 7+0.7*i;
% end
% radiomap_index = 65:(65+8);             % n1 2f east

%% KI 1F
parentsFolders = getNameFolds('rawdata/ki*');
b = cell(length(parentsFolders),1);

for i = 1:length(parentsFolders)
     a = getNameFolds(sprintf('rawdata/%s',parentsFolders{i}));
     b{i} = cellfun(@(x) sprintf('rawdata/%s/%s',parentsFolders{i},x), a,...
         'UniformOutput',false);
%      b = what(sprintf('rawdata/%s/%s',parentsFolders{i},a{1}));
end
b = cat(1,b{:});       % all path info from rawdata (radiomap survey lines)

path_coord = cell(1,length(b));
for i = 1:3
    path_coord{i}.x = [3,44.4];
    path_coord{i}.y = 7.5+0.6*(i-1);
end
for i = 4:4+13
    path_coord{i}.x = [3,9.8];
    path_coord{i}.y = 7.5+0.6*(i-1);
end
for i = 18:18+1
    path_coord{i}.x = 3+0.9*(i-18);
    path_coord{i}.y = [17.2,29.8];
end
for i = 20:20+1
    path_coord{i}.x = [3.9,7.6];
    path_coord{i}.y = 30.4+0.6*(i-20);
end
for i = 22:22+4
    path_coord{i}.x = [5.7,10.2+1.8*(i-22)];
    path_coord{i}.y = 31.5+0.6*(i-22);
end
for i = 27:27+1
    path_coord{i}.x = [5.7,21.9];
    path_coord{i}.y = 34.5+0.6*(i-27);
end
for i = 29:29+3
    path_coord{i}.x = [16.5,21.9];
    path_coord{i}.y = 35.7+0.6*(i-29);
end
for i = 33:33+11
    path_coord{i}.x = [31.8,44.4];
    path_coord{i}.y = 9.3+0.6*(i-33);
end
for i = 45:45+30
    path_coord{i}.x = [31.8,50.7];
    path_coord{i}.y = 16.5+0.6*(i-45);
end
for i = 76:76+18
    path_coord{i}.x = [21.9,50.7];
    path_coord{i}.y = 35.1+0.6*(i-76);
end
%%
% for i = 1:length(radiomap_index)
for i = 1:length(b)
%     target_rawdata_paths = getNameFolds('rawdata/ki*');
%     rawdata = load_rawdata(fullfile('rawdata',target_rawdata_paths{radiomap_index(i)}));
    
%     rawdata = load_rawdata(fullfile('rawdata',target_rawdata_paths{i+t_input_idx}));
    rawdata = load_rawdata(b{i});

    %% resample for synchronize
    rate = 2e-2;
    processed_data = resample_rawdata(rawdata,rate);

    %% find step point (step) and time labeling
    % threshold should be tuned experimentally to match a person's level 
    minPeakHeight = std(processed_data.acc_norm);       
    [~,locs] = findpeaks(processed_data.acc_norm,'MinPeakDistance',...
        .3/rate,'MinPeakHeight',minPeakHeight);   % .3s 이내의 피크는 무시 (가정: 발걸음이 .3초 내에는 2번이뤄지지 않음)\

    if length(path_coord{i}.x) == 2     % case of horizontal path
        gt_x = linspace(path_coord{i}.x(1),path_coord{i}.x(2),length(locs));
        gt_y = path_coord{i}.y(1)*ones(1,length(locs));
    elseif length(path_coord{i}.y) == 2 % case of vertical path
        gt_x = path_coord{i}.x(1)*ones(1,length(locs));
        gt_y = linspace(path_coord{i}.y(1),path_coord{i}.y(2),length(locs));
    end
    
    if ~exist('map','var')
        map = [gt_x',gt_y',processed_data.Magnetometer(locs,:)];
    else
        map = [map;gt_x',gt_y',processed_data.Magnetometer(locs,:)];
    end
    %% v2
%     d_locs = diff(locs);
% %     r_gt_x = [];
% %     r_gt_y = [];
%     remap = [];
%     for j = 1:length(d_locs)
%         r_gt_x = linspace(gt_x(j),gt_x(j+1),d_locs(j)+1);
%         r_gt_y = linspace(gt_y(j),gt_y(j+1),d_locs(j)+1);
%         r_gt_x(end) = [];
%         r_gt_y(end) = [];
%         remaped_mag = processed_data.Magnetometer(locs(j):(locs(j)+d_locs(j)-1),:);
%         if isempty(remap)
%             remap = [r_gt_x',r_gt_y',remaped_mag];
%         else
%             remap = [remap; r_gt_x',r_gt_y',remaped_mag];
%         end
%     end
%     
%     % processed_data.Magnetometer(locs,:)
%     if ~exist('map','var')
% %         map = [gt_x',gt_y',processed_data.Magnetometer(locs,:)];
%         map = remap;
%     else
% %         map = [map;gt_x',gt_y',processed_data.Magnetometer(locs,:)];
%         map = [map; remap];
%     end
end
mat_filename = fullfile('mats',['magmap-ki-1f-',num2str(0.6),'p.mat']);
if ~exist(mat_filename,'file')
    save(mat_filename,'map')
end
