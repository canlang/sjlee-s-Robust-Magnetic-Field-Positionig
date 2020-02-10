function map = magmap_construction(mypath,site_name,d)
% d: interval
% addpath(genpath(path));
if exist(fullfile(mypath,['magmap-',site_name,num2str(d),'a.mat']), 'file') == 2
    load(fullfile(mypath,['magmap-',site_name,num2str(d),'a.mat']), 'map')    
else
    
    switch site_name
        case 'N1-7F'
            data1 = readtable('batch.csv');    
            lM = [data1.magnet_x,data1.magnet_y,data1.magnet_z];
            % INTERPOLATION
            x = data1.x;
            y = data1.y;
        case 'KI-1F'
            data = load(sprintf('mats/magmap-%s-0.6p.mat',site_name));
            lm.x = data.map(:,1);lm.y = data.map(:,2);
            lM = data.map(:,3:5);
            x = lm.x;
            y = lm.y;
    end
    newlM = [];
    for i=1:3
        z = lM(:,i);
        XI = min(x):d:max(x);
        YI = min(y):d:max(y);
        [X,Y] = meshgrid(XI,YI);
        shp = alphaShape(x,y);
        
        disp(shp.Alpha)
        shp.Alpha = .6;
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
    save(fullfile(mypath,['magmap-',site_name,num2str(d),'a.mat']),'map')
end

