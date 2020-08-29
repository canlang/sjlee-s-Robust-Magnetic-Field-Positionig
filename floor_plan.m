clearvars;close all

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
    path_coord{i}.x = 3+0.9*(i-18)+0.4;
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

site_name = 'KI-1F';
% site_name = "N1-7F";
A = imread(sprintf('map/%s-2.png',site_name),'BackgroundColor',[1 1 1]);
switch site_name    
    case 'KI-1F'
        scale = 0.05;
        xWorldLimits = [0 size(A,2)*scale];
        yWorldLimits = [0 size(A,1)*scale];
    case 'N1-7F'        
        xWorldLimits = [-1 1650/20];
        yWorldLimits = [-1 660/20];
end

%%
RA = imref2d(size(A),xWorldLimits,yWorldLimits);
imshow(flipud(A),RA);
axis xy
hold on
for i=1:length(path_coord)
    if length(path_coord{i}.x) == 2
        plot(path_coord{i}.x,path_coord{i}.y*ones(2,1),'b-')
    else
        plot(path_coord{i}.x*ones(2,1),path_coord{i}.y,'b-')
    end
end
hold off
% xlim([3.0000-5   50.5000+10])
% ylim([7.5000-5   45.5000+5])
legend('data collecting path')    
sdf(gcf,'sj')