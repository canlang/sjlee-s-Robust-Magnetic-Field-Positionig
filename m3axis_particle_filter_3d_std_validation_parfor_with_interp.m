clear; close all;
data1 = readtable('batch.csv');
data2 = readtable('20171124 MagCoord3axisData.csv');
lM = [data1.magnet_x,data1.magnet_y,data1.magnet_z];

% INTERPOLATION
interp_interval = .8;
x = data1.x;
y = data1.y;
newlM = [];
for i=1:3
    z = lM(:,i);
    XI = min(x):interp_interval:max(x);
    YI = min(y):interp_interval:max(y);
    [X,Y] = meshgrid(XI,YI);
    shp = alphaShape(x,y);
    in = inShape(shp,X,Y);
    xg = X(in);
    yg = Y(in);
    zg = griddata(x,y,z,xg,yg);         % 2. griddata() : INTERPOLATION
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
lM = newlM(:,3:end);
%%
% nParticleCandidate = 500:500:3000;
% nRepeat = 500;
% nParticleCandidate = 2000:100:3000;
% nRepeat = 500;
nParticleCandidate = 1000:1000:3000;
% % nParticleCandidate = [1000 2000];
nRepeat = 100;

% errMat = zeros(length(nParticleCandidate),nRepeat);
% stdMat = zeros(length(nParticleCandidate),nRepeat);
convIndexes = zeros(length(nParticleCandidate),nRepeat);    
errMat = cell(length(nParticleCandidate),nRepeat);
%%
for k = 1:length(nParticleCandidate)
    parfor j = 1:nRepeat
        % ------------------------
        % initialize particle
        n = nParticleCandidate(k);
        
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
            [phy_dist,I] = pdist2([data1.x,data1.y],[ps_x,ps_y],'euclidean','Smallest',1);
            % 2. calculate magnetic distance
            R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]),ps_mag_heading,'UniformOutput',false);
            rotatedMag = cell2mat(cellfun(@(x)((x*tM(i,:)')'),R,'UniformOutput',false));
            % TODO: may be more optimizable (DONE?maybe)
        %     mag_dist = diag(pdist2(rotatedMag,lM(I,:),'euclidean'));
            mag_dist = sqrt(sum((rotatedMag-lM(I,:)).^2,2));
        %     mag_dist = bsxfun(@(x,y) pdist([x;y]), rotatedMag,lM(I,:));

            if all(mag_dist) == 0
                break
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
        convIndexes(k,j) = cIdx;
        if cIdx<(length(data2.coord_x)/2)
            errMat{k,j} = est_err(cIdx:end);
        end
        
%         errMat{k,j} = est_err;
        
        
%         errMat(k,j) = mean(mean(err(63:end,:)));
%         stdMat(k,j) = std(std(err(63:end,:)));
    end
end
%%
close all
figure
% axes('ColorOrder',brewermap(3,'RdYlGn'),'NextPlot','replacechildren')
subplot(1,3,1)
boxplot(convIndexes',nParticleCandidate)
xlabel('Number of particles N')
ylabel('Distance to convergence (m)')
set(gca,'XTickLabelRotation',45)

subplot(1, 3, [2 3])

% colormap(brewermap([],''))
mean_err = zeros(length(nParticleCandidate),1);

hold on
linS = {'--','-.','-',':','--'};
for i = 1:length(nParticleCandidate)
    h = cdfplot(cell2mat(errMat(i,:)'));
    set(h,'lineStyle',linS{i})
%     set(h,'lineStyle',linS{i},'linewidth',4)
    mean_err(i) = mean(cell2mat(errMat(i,:)'));
end
txt = strcat({'N='},int2str(nParticleCandidate'),num2str(mean_err,'\t ,(Mean=%.2f)'));
legend(txt,'location','best')
xlabel('Error distance (m)')
set(0,'DefaultAxesColorOrder',brewermap(3,'Set1'))

set(gcf,'units','points','position',[200,500,1100,500])
% sdf(gcf,'sj3')
tightfig(gcf)

save(sprintf('est-result/n1-7f-parfor-%s.mat',num2str(interp_interval)),'convIndexes','errMat','nParticleCandidate')
return
%%
% plot(nParticleCandidate,mean(errMat,2),'x-')

boxPlot = @iosr.statistics.boxPlot;
% boxPlot(errMat','style','hierarchy','groupLabels',num2cell(string(nParticleCandidate)))
% boxPlot(errMat','style','hierarchy','groupLabels',nParticleCandidate)
figure('color','w');
h2 = boxPlot(nParticleCandidate,errMat', ...
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
notBoxPlot(errMat(5:end,:)',nParticleCandidate(5:end),'style','line','markMedian',true,'interval','SEM')