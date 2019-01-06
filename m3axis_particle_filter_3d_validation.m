clear; close all;
data1 = readtable('batch.csv');
data2 = readtable('20171124 MagCoord3axisData.csv');
lM = [data1.magnet_x,data1.magnet_y,data1.magnet_z];

%%
N = 5;
for j = 1:N
    % initialize particle
    n = 1500;
    % 1. only road
    rand_idx = randi(length(data1.x),n,1);
    ps.x = data1.x(rand_idx);
    ps.y = data1.y(rand_idx);

    % 2. all area
    % ps.x = random('Uniform', min(data1.x),max(data1.x),n,1);
    % ps.y = random('Uniform', min(data1.y),max(data1.y),n,1);

    % 3. initial area
    % ps.x = data2.coord_x(1)+random('normal',0,5,n,1);
    % ps.y = data2.coord_y(1)+random('normal',0,5,n,1);
    % ps.y = 15+random('normal',0,5,n,1);

    ps.mag_heading = random('Uniform', 0,2*pi,n,1);
    ps.phy_heading = random('Uniform', 0,2*pi,n,1);
    ps.prob = ones(n,1)*(1/n);
    % ps.stlng = ones(n,1) + random('Uniform', -.1,.1,n,1);


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
    %     ps.x = bsxfun(@(x,y) x + cos(y),ps.x,ps.phy_heading);
    %     ps.y = bsxfun(@(x,y) x + sin(y),ps.y,ps.phy_heading);
    %     ps.x = ps.x + cos(ps.mag_heading-pi/2);     % when heading shifted input date
    %     ps.y = ps.y + sin(ps.mag_heading-pi/2);
        ps.x = ps.x + cos(ps.mag_heading);
        ps.y = ps.y + sin(ps.mag_heading);
    %     ps.x = ps.x + ps.stlng.*cos(ps.heading);
    %     ps.y = ps.y + ps.stlng.*sin(ps.heading);

    % ================ update

        % 1. find (geo-locational) nearest learning data
        [phy_dist,I] = pdist2([data1.x,data1.y],[ps.x,ps.y],'euclidean','Smallest',1);
        % 2. calculate magnetic distance
        R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]),-ps.mag_heading,'UniformOutput',false);
        rotatedMag = cell2mat(cellfun(@(x)((x*tM(i,:)')'),R,'UniformOutput',false));
        % TODO: may be more optimizable (DONE?maybe)
    %     mag_dist = diag(pdist2(rotatedMag,lM(I,:),'euclidean'));
        mag_dist = sqrt(sum((rotatedMag-lM(I,:)).^2,2));
    %     mag_dist = bsxfun(@(x,y) pdist([x;y]), rotatedMag,lM(I,:));

        if all(mag_dist) == 0
            break
        end
        ps.prob = 1./(mag_dist);
        ps.prob(phy_dist>1) = 0;
        if sum(ps.prob) == 0
            rand_idx = randi(length(data1.x),n,1);
            ps.x = data1.x(rand_idx);
            ps.y = data1.y(rand_idx);
            ps.mag_heading = random('Uniform', 0,2*pi,n,1);
            ps.phy_heading = random('Uniform', 0,2*pi,n,1);
            ps.prob = ones(n,1)*(1/n);
        else
            ps.prob = ps.prob./sum(ps.prob);
        end

        % ================ resample
        resample_idx = randsample(1:n,n,true,ps.prob);
        phy_move_noise_range = 1;
        ps.x = ps.x(resample_idx) + phy_move_noise_range*rand(n,1) - phy_move_noise_range/2;
        ps.y = ps.y(resample_idx) + phy_move_noise_range*rand(n,1) - phy_move_noise_range/2;
        ps.mag_heading = ps.mag_heading(resample_idx)+random('normal',0,.001,n,1);
    %     ps.phy_heading = ps.phy_heading(resample_idx)+random('normal',0,.001,n,1);

    %     ps.heading = ps.heading(resample_idx)+random('Uniform', -pi/10,pi/10,n,1);
    %     ps.stlng = ps.stlng(resample_idx) + random('normal',0,1,n,1);
    %     ps.prob = init_prob;

    
        % Particle's summary
    %     est(i,:) = [mean(ps.x),mean(ps.y),mean(angdiff(ps.mag_heading,0))];
        est(i,:) = [mean(ps.x),mean(ps.y),mean(abs(angdiff(ps.mag_heading,0)))];        %x,y,heading
        
        % calc Euclidiean distance for error validation
        err(i,:) = pdist2([data2.coord_x(i),data2.coord_y(i)],[ps.x,ps.y]);
    
    %     break
    %     pause(.1)

    end
    delete(subplot(N,1,j))
    subplot(N,1,j)
    convergence_thr = pdist2([data2.coord_x(end),data2.coord_y(end)],est(end,1:2));
    if convergence_thr < 1
        [lineh, bandsh] = fanChart(1:size(err,1),err, 'mean', 10:10:90, ...
            'alpha', .2, 'colormap', {'shadesOfColor', [0 0 .8]});
    else
        [lineh, bandsh] = fanChart(1:size(err,1),err, 'mean', 10:10:90, ...
            'alpha', .2, 'colormap', {'shadesOfRed'});
    end
    txt = strcat({'Pct'}, cellstr(int2str((20:20:80)')));
    legend([bandsh;lineh], [txt;{'Mean'}])
    xlim([0,63])
%     errs{j} = err;

%     clickableLegend([lineh;bandsh], [{'Mean'};txt])
end

%%
% samexaxis('abc','xmt','on','ytac','join','yld',1)
suplabel('Location index (time flow)');
suplabel('Error distance (m)','y');

% samexaxis('xmt','off','ytac','yld',1,'box','on','join')
set(gcf,'units','points','position',[1200,500,800,800])
sdf(gcf,'sj2')

print -clipboard -dbitmap