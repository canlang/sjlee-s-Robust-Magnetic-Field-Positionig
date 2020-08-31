% clear;
close all;
% purpose: just evaluation and plot, temporal
% copied from: m3axis_particle_filter_3d_std_validation_new_inputform_v2

site_name = 'N1-7F';
atti = [85 45 0 20];
repeat = 50;
errs = {repeat,length(atti)};

for j=1:length(atti)
    parfor k=1:repeat
        t_input_idx = 90+(j-1);   

        intp_intv = .6;
        map = loadMagneticMap('mats',site_name,intp_intv);

        lm_x = map(:,1);lm_y = map(:,2);
        lM = map(:,3:5);


        target_rawdata_paths = getNameFolds('rawdata');
        fprintf('loading %s data...\n',target_rawdata_paths{t_input_idx});
        rawdata = load_rawdata(fullfile('rawdata',target_rawdata_paths{t_input_idx}));
        %% resample for synchronize
        rate = 2e-2;
        % TODO: refactoring
        processed_data = resample_rawdata(rawdata,rate);          
        % differ: gryo radian, time column
        % processed_data = resample_rawdata2(rawdata,rate);
        % processed_data.Gyroscope = processed_data.Gyroscope*(180/pi);

        %% find step point (step) and time labeling
        % threshold should be tuned experimentally to match a person's level 
        minPeakHeight = std(processed_data.acc_norm);       
        [pks,locs] = findpeaks(processed_data.acc_norm,'MinPeakDistance',...
            .3/rate,'MinPeakHeight',minPeakHeight);   % .3s 이내의 피크는 무시 (가정: 발걸음이 .3초 내에는 2번이뤄지지 않음)\

        %%
        addpath(genpath('madgwick_algorithm_matlab'));
        AHRS = MadgwickAHRS('SamplePeriod', rate, 'Beta', 0.1); % sample rate: 2e-2

        quaternion = zeros(length(processed_data.Time), 4);
        for t = 1:length(processed_data.Time)
            AHRS.Update(processed_data.Gyroscope(t,:) * (pi/180), ...
                processed_data.Accelerometer(t,:), ...
                processed_data.Accelerometer(t,:));	% gyroscope units must be radians
            quaternion(t, :) = AHRS.Quaternion;
        end

        % euler : use conjugate for sensor frame relative to Earth and convert to degrees.
        euler = quatern2euler(quaternConj(quaternion(locs,:)));	
        % euler = quatern2euler(quaternConj(quaternion)) * (180/pi);
        d_heading = [0;diff(euler(:,3))];

        rotMat = quatern2rotMat(quaternion(locs,:));
        % rotMat = quatern2rotMat(quaternConj(quaternion(locs,:)));

        std_euler = stdfilt(unwrap(euler(:,3)));
        % std_euler = stdfilt(deg2rad(unwrap(euler(:,3))));
        % plot(processed_data.Time, std_euler)

        %% inbound & outbound
        % layout = loadjson('map/N1-7F-HiRes2.json');
        layout = jsondecode(fileread(sprintf('map/%s.json',site_name)));

        x = layout.in(:,1);
        y = layout.in(:,2);
        shp = polyshape(x,y);
        for i = 1:length(layout.out)
            ox = layout.out{i}(:,1);
            oy = layout.out{i}(:,2);
            shp = subtract(shp,polyshape(ox,oy));
        end

        A = imread('map/N1-7F.png','BackgroundColor',[1 1 1]);

        xWorldLimits = [-1 1650/20];
        yWorldLimits = [-1 660/20];
        RA = imref2d(size(A),xWorldLimits,yWorldLimits);
        imshow(flipud(A),RA);
        axis xy;


        % draw learning data
        hold on
        plot(lm_x,lm_y,'.','MarkerSize', 10)
        % plot(lM(:,1),lM(:,2),'.','MarkerSize', 10)
        % for save eps
        legend('reference point')
        sdf(gcf,'sj2')
        % print -depsc2 env_setting.eps

        % axis equal
        % xlim([8 83])
        % ylim([0 30])
        % set(gcf,'units','points','position',[700,500,1500,700])
        plot(shp,'FaceAlpha',.5,'EdgeColor','r')


        %%
        % initialize particle
        n = 3000;
        % 1. only road
        rand_idx = randi(length(lm_x),n,1);
        ps_x = lm_x(rand_idx);
        ps_y = lm_y(rand_idx);

        % 2. all area
        % ps_x = random('Uniform', min(lm_x),max(lm_x),n,1);
        % ps_y = random('Uniform', min(lm_y),max(lm_y),n,1);

        % 3. initial area
        % % ps_x = data2.coord_x(1)+random('normal',0,5,n,1);
        % % ps_y = data2.coord_y(1)+random('normal',0,5,n,1);
        % ps_x = 7+random('normal',0,5,n,1);
        % ps_y = 12+random('normal',0,5,n,1);

        ps_sl = ones(n,1)*.7;
        ps_mag_heading = random('Uniform', 0,2*pi,n,1);
        % ps_phy_heading = random('Uniform', 0,2*pi,n,1);
        ps_prob = ones(n,1)*(1/n);
        % ps_stlng = ones(n,1) + random('Uniform', -.1,.1,n,1);

        % draw particle
        hold on 
        % h_ps = plot(ps_x,ps_y,'.','MarkerSize',8);
        h_ps = scatter(ps_x,ps_y,20,'c','filled','MarkerFaceAlpha',.2);
        h_pm = plot(mean(ps_x),mean(ps_y),'ms');
        % h_gt = plot(data2.coord_x(1),data2.coord_y(1),'s','MarkerSize',10,...
        %     'MarkerEdgeColor','b',...
        %     'MarkerFaceColor',[0.5,0.5,0.5]);
        hold off


        % test data matrix
        % tM = [data2.mag_x,data2.mag_y,data2.mag_z];
        % rotate test data 
        % x = -pi/2;
        x = 0;
        R = [cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1];

        % tM = (R*tM')';
        % tM = (R*tM')'-25;
        % 25 offset meaning TestData shift: because magnetometer's calibration not matched
        tM = processed_data.Magnetometer(locs,:);

        % test result matrix
        est = zeros(length(tM),3);
        AD = zeros(length(tM),n);
        eta = zeros(length(tM),n);      % heading diff 
        %%
        % yaw_offset = pi/2;       % pi/2
        yaw_offset = 0;       % pi/2

        for i = 1:length(tM)
            % ================ PREDICTION
        %     ps_x = bsxfun(@(x,y) x + cos(y),ps_x,ps_phy_heading);
        %     ps_y = bsxfun(@(x,y) x + sin(y),ps_y,ps_phy_heading);
        %     ps_x = ps_x + cos(ps_mag_heading-pi/2);     % when heading shifted input date
        %     ps_y = ps_y + sin(ps_mag_heading-pi/2);
            hN = .03;      % heading noise value
            if i>1
        %         Halpha = pi-est(i-1,3);
                Halpha = 0;
                ps_mag_heading = ps_mag_heading+random('normal',Halpha,hN,n,1);     % candidate, .08
            else
        %     ps_mag_heading = ps_mag_heading+euler(i,3);
                ps_mag_heading = ps_mag_heading+random('normal',0,hN,n,1);     % candidate, .08
            end

            ps_sl = .7 + random('normal',0,.5,n,1);

            mu = [0,0];        

            ps_x = ps_x + cos(ps_mag_heading+euler(i,3)).*ps_sl+ random('Uniform',-.1,.1,n,1);
            ps_y = ps_y + sin(ps_mag_heading+euler(i,3)).*ps_sl+ random('Uniform',-.1,.1,n,1);

            % ================ UPDATE
            % 1. find (geo-locational) nearest learning data
            [phy_dist,I] = pdist2([lm_x,lm_y],[ps_x,ps_y],'euclidean','Smallest',1);
            % 2. calculate Rotated magnetic field data and magnetic distance

            % ILoA, line 256, 260, 276
            rotatedMag2 = getHeadingRotatedVector(ps_mag_heading, tM(i,:), rotMat(:,:,i));
            rotatedMag = rotatedMag2;
            mag_dist = sqrt(sum((rotatedMag-lM(I,:)).^2,2));

            % COSINE
        %     mag_dist = diag(pdist2(rotatedMag,lM(I,:),'minkowski',3));


            % MALOC
        %     observed = [vecnorm(lM(I,1:2),2,2), lM(I,3)];
        %     measured = [vecnorm(tM(i,1:2),2), tM(i,3)];
        %     mag_dist2 = pdist2(observed, measured,'mahalanobis',nearestSPD(nancov(observed)));
        %     mag_dist = mag_dist2';

            if ~all(mag_dist)
                break
            end
            ps_prob = 1./(mag_dist);
            in = isinterior(shp,ps_x,ps_y);
        %     in = inShape(shp,[ps_x,ps_y]);
            ps_prob(~in) = 0;
            ps_prob(phy_dist'>3) = 0;

            if sum(ps_prob) == 0
                rand_idx = randi(length(lm_x),n,1);
                ps_x = lm_x(rand_idx);
                ps_y = lm_y(rand_idx);
                ps_mag_heading = random('Uniform', 0,2*pi,n,1);
        %         ps_phy_heading = random('Uniform', 0,2*pi,n,1);
                ps_prob = ones(n,1)*(1/n);
            else
                ps_prob = ps_prob./sum(ps_prob);
            end

            % ================ RESAMPLE
            resample_idx = randsample(1:n,n,true,ps_prob);
            ps_x = ps_x(resample_idx);
            ps_y = ps_y(resample_idx);
            ps_mag_heading = ps_mag_heading(resample_idx);
            ps_sl = ps_sl(resample_idx);

            set(h_ps,'XData',ps_x,'YData',ps_y)                             % ps result
            set(h_pm,'XData',mean(ps_x),'YData',mean(ps_y));
        %     set(h_gt,'XData',data2.coord_x(i),'YData',data2.coord_y(i))     % ground truth
            drawnow

            est(i,:) = [mean(ps_x),mean(ps_y),circ_mean(ps_mag_heading)]; 
            eta(i,:) = wrapTo2Pi(ps_mag_heading'); 

        %     err(i,:) = pdist2([data2.coord_x(i),data2.coord_y(i)],[ps_x,ps_y]);
            AD(i,:) = pdist2([mean(ps_x),mean(ps_y)],[ps_x,ps_y]);
        end


        std_AD = std(AD,0,2);
        converge_idx = find(std_AD <= 2,1);
        % heading_err = pdist2(euler(:,3),est(:,3));
        errs{k,j} = abs(euler(converge_idx:end,3)+est(converge_idx:end,3));
    end
end
%%
close all
figure
x1 = vertcat(errs{1})/pi*180; x2 = vertcat(errs{2})/pi*180; x3 = vertcat(errs{3})/pi*180;
x4 = vertcat(errs{4})/pi*180;
x = [x3;x4;x2;x1];
g = [ones(size(x3)); 2*ones(size(x4)); 3*ones(size(x2)); 4*ones(size(x1))];
boxplot(x,g)
fprintf('0deg:%.2f, 20deg:%.2f, 45deg:%.2f, 90deg:%.2f\n',mean(x3),mean(x4),mean(x2),mean(x1));

figure
hold on
for i=1:2
    cdfplot(vertcat(errs{:,i})/pi*180)
end
cdfplot(vertcat(errs{:,4})/pi*180)
cdfplot(vertcat(errs{:,3})/pi*180)
legend('85','45','20','0')