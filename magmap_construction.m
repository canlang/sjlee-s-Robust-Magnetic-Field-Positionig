function map = magmap_construction(path,d)
% d: interval
% addpath(genpath(path));
if exist(fullfile(path,['magmap',num2str(d),'.mat']), 'file') == 2
    load(fullfile(path,['magmap',num2str(d),'.mat']), 'map')    
else
    data1 = readtable('batch.csv');
    lM = [data1.magnet_x,data1.magnet_y,data1.magnet_z];
    % INTERPOLATION
    x = data1.x;
    y = data1.y;
    newlM = [];
    for i=1:3
        z = lM(:,i);
        XI = min(x):d:max(x);
        YI = min(y):d:max(y);
        [X,Y] = meshgrid(XI,YI);
        shp = alphaShape(x,y);
        in = inShape(shp,X,Y);
        xg = X(in);
        yg = Y(in);
        zg = griddata(x,y,z,xg,yg,'nearest');         % 2. griddata() : INTERPOLATION
        if isempty(newlM)
            newlM = [xg,yg,zg];
        else
            newlM = [newlM,zg];
            if ~isequal(newlM(:,1:2), [xg,yg])
                disp 'error'
            end
        end
    end
    data1 = array2table(newlM(:,1:2), 'VariableNames',{'x','y'});
    map = [data1.x, data1.y, newlM(:,3:end)];
    save(fullfile(path,['magmap',num2str(d),'.mat']),'map')
end

