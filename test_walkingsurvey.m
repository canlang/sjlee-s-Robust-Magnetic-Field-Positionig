clear all

path_coord{1}.x = [9,70.7];
path_coord{1}.y = 21;

path_coord{2}.x = [12,80];
path_coord{2}.y = 8;

path_coord{3}.x = 22;
path_coord{3}.y = [21,8];

path_coord{4}.x = 38.3;
path_coord{4}.y = [21,8];

path_coord{5}.x = 54;
path_coord{5}.y = [21,8];

path_coord{6}.x = 71;
path_coord{6}.y = [21,8];

path_coord{7}.x = [9,70.7];
path_coord{7}.y = 21.6;

path_coord{8}.x = [9,70.7];
path_coord{8}.y = 20.4;

path_coord{9}.x = [12,80];
path_coord{9}.y = 8.6;

path_coord{10}.x = [12,80];
path_coord{10}.y = 7.4;
%%

% t_input_idx = 43;       % start index for radio-map rawdata.
% radiomap_index = [44:(44+5),57:(57+3)];     % east-direction
radiomap_index = [64];
for i = 1:length(radiomap_index)
    target_rawdata_paths = getNameFolds('rawdata');
    rawdata = load_rawdata(fullfile('rawdata',target_rawdata_paths{radiomap_index(i)}));
%     rawdata = load_rawdata(fullfile('rawdata',target_rawdata_paths{i+t_input_idx}));


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
save(fullfile('mats',['magmap',num2str(0.6),'.mat']),'map')

function nameFolds = getNameFolds(pathFolder)
        d = dir(pathFolder);
        isub = [d(:).isdir]; %# returns logical vector
        nameFolds = {d(isub).name}';
        nameFolds(ismember(nameFolds,{'.','..'})) = [];
        % nameFolds(~cellfun(@isempty,regexp(nameFolds, '\d{6}'))) = [];
end
    