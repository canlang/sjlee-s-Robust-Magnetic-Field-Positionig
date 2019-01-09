% clear; close all;
% learning data
data1 = readtable('batch.csv');
lM = [data1.magnet_x,data1.magnet_y,data1.magnet_z];
% testing data
% data2 = readtable('20171124 MagCoord3axisData.csv');
target_rawdata_paths = getNameFolds('rawdata');
rawdata = load_rawdata(fullfile('rawdata',target_rawdata_paths{1}));
trace_info = regexp(target_rawdata_paths{j},'_','split');

%%
% nParticleCandidate = 500:500:3000;
% nRepeat = 500;
% nParticleCandidate = 2000:100:3000;
% nRepeat = 500;
% nParticleCandidate = 500:250:3000;
% nRepeat = 100;
nParticleCandidate = 500:250:2000;
nRepeat = 5;

errMat = zeros(length(nParticleCandidate),nRepeat);
stdMat = zeros(length(nParticleCandidate),nRepeat);
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
            R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]),-ps_mag_heading,'UniformOutput',false);
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
            phy_move_noise_range = 1;
            ps_x = ps_x(resample_idx) + phy_move_noise_range*rand(n,1) - phy_move_noise_range/2;
            ps_y = ps_y(resample_idx) + phy_move_noise_range*rand(n,1) - phy_move_noise_range/2;
            ps_mag_heading = ps_mag_heading(resample_idx)+random('normal',0,.001,n,1);
        %     ps_phy_heading = ps_phy_heading(resample_idx)+random('normal',0,.001,n,1);

        %     ps_heading = ps_heading(resample_idx)+random('Uniform', -pi/10,pi/10,n,1);
        %     ps_stlng = ps_stlng(resample_idx) + random('normal',0,1,n,1);
        %     ps_prob = init_prob;

        %     est(i,:) = [mean(ps_x),mean(ps_y),mean(angdiff(ps_mag_heading,0))];
        
            est(i,:) = [mean(ps_x),mean(ps_y),mean(abs(angdiff(ps_mag_heading,0)))];        %x,y,heading
            err(i,:) = pdist2([data2.coord_x(i),data2.coord_y(i)],[ps_x,ps_y]);
        end
        
%         errMat(k,j) = mean(err(50,:));
%         stdMat(k,j) = std(err(50,:));
        
%         errMat(k,j) = mean(mean(err));
%         stdMat(k,j) = std(std(err));
        
        errMat(k,j) = mean(mean(err(63:end,:)));
        stdMat(k,j) = std(std(err(63:end,:)));
    end
end
%%
boxplot(errMat(5:end,:)',nParticleCandidate(5:end))

% boxplot(errMat3',nParticleCandidate3)
% set(gca,'XTickLabelRotation',45)

xlabel('N (number of particles)')
ylabel('Error distance (m)')
% sdf(gcf,'sj2')
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

%% local functions -----------------------------------------------------------------
% 1st local function
function nameFolds = getNameFolds(pathFolder)
d = dir(pathFolder);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];
% nameFolds(~cellfun(@isempty,regexp(nameFolds, '\d{6}'))) = [];
end

% 2nd local function
function data = resample_rawdata(rawdata,rate)
% Var1 : 3 vectors of Accelerometer, Var2: norm vector of Acc.
T_acc = timetable(seconds(rawdata.acc(:,2)/1e9),rawdata.acc(:,3:5),rawdata.acc_norm,...
    'VariableNames',{'acc','acc_norm'});
% Var1 : 3 vectors of Gyroscope
T_gyr = timetable(seconds(rawdata.gyr(:,2)/1e9),rawdata.gyr(:,3:5),...
    'VariableNames',{'gyr'});
T_acc = sortrows(T_acc);
T_gyr = sortrows(T_gyr);

TT = synchronize(T_acc,T_gyr,'regular','linear','TimeStep',seconds(rate));

data.Accelerometer = TT.acc;
data.Gyroscope = TT.gyr*180/pi;
data.acc_norm = TT.acc_norm;

data.Time = seconds(TT.Time(:));
% data.Time = seconds(TT.Time(:)-(TT.Time(1)));
data.Rate = median(diff(data.Time)); % cal sample rate
end