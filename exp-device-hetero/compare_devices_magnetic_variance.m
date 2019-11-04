clearvars;close all;
target_rawdata_paths = getNameFolds('./');
data1 = csvread(fullfile('./',target_rawdata_paths{1},'MagneticField.csv'),1,0);
data2 = csvread(fullfile('./',target_rawdata_paths{2},'MagneticField.csv'),1,0);

temp_data1 = vecnorm(data1(:,3:5),2,2);
temp_data2 = vecnorm(data2(:,3:5),2,2);

plot3(data1(:,3),data1(:,4),data1(:,5),'o')
hold on
plot3(data2(:,3),data2(:,4),data2(:,5),'o')
% axis equal 
% pbaspect([1 1 1])
function nameFolds = getNameFolds(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
% nameFolds(~cellfun(@isempty,regexp(nameFolds, '\d{6}'))) = [];
end
