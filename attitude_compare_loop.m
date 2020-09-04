% clear;
close all;
% purpose: just evaluation and plot, temporal
% copied from: m3axis_particle_filter_3d_std_validation_new_inputform_v2

site_name = 'N1-7F';
atti = [85 45 0 20];
repeat = 10;
errs = {repeat,length(atti)};

for j=1:length(atti)
    for k=1:repeat
        t_input_idx = 90+(j-1);   

%     %     gt_filename = sprintf('mats/magmap-n1-7f-step-wpNCE-uturn%d copy.mat',atti(t_input_idx-86));
%         gt_filename = sprintf('mats/magmap-n1-7f-step-wpNCE%d.mat',atti(j));
%         load(gt_filename,'map');
%         err = zeros(length(map),1);
%         gt.x = map(:,1);gt.y = map(:,2);

        intp_intv = .6;
        map = loadMagneticMap('mats',site_name,intp_intv);

        lm.x = map(:,1);lm.y = map(:,2);
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
        plot(lm.x,lm.y,'.','MarkerSize', 10)
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
        n = 2000;
        % 1. only road
        rand_idx = randi(length(lm.x),n,1);
        ps.x = lm.x(rand_idx);
        ps.y = lm.y(rand_idx);

        % 2. all area
        % ps.x = random('Uniform', min(lm.x),max(lm.x),n,1);
        % ps.y = random('Uniform', min(lm.y),max(lm.y),n,1);

        % 3. initial area
        % % ps.x = data2.coord_x(1)+random('normal',0,5,n,1);
        % % ps.y = data2.coord_y(1)+random('normal',0,5,n,1);
        % ps.x = 7+random('normal',0,5,n,1);
        % ps.y = 12+random('normal',0,5,n,1);

        ps.sl = ones(n,1)*.7;
        ps.mag_heading = random('Uniform', 0,2*pi,n,1);
        % ps.phy_heading = random('Uniform', 0,2*pi,n,1);
        ps.prob = ones(n,1)*(1/n);
        % ps.stlng = ones(n,1) + random('Uniform', -.1,.1,n,1);

        % draw particle
        hold on 
        % h_ps = plot(ps.x,ps.y,'.','MarkerSize',8);
        h_ps = scatter(ps.x,ps.y,20,'c','filled','MarkerFaceAlpha',.2);
        h_pm = plot(mean(ps.x),mean(ps.y),'ms');
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
        switch t_input_idx
            case 43
                video_filename = 'm3axis_2d-space_pf_real-time-input_loop-traj';        %input index = 43;
            otherwise
        %         video_filename = 'm3axis_2d_pf_nonshiftedinput_newinput_rotated_6';%nonshifted input??
                video_filename = target_rawdata_paths{t_input_idx};
        end
        switch video_flag
            case 1
                v = VideoWriter(strcat('vids/',video_filename,'.mp4'),'MPEG-4');
                v.FrameRate = 10;
                v.Quality = 100;
                open(v);
                frame = getframe(gcf);
                writeVideo(v,frame);
        end

        % yaw_offset = pi/2;       % pi/2
        yaw_offset = 0;       % pi/2

        for i = 1:length(tM)
            % ================ PREDICTION
        %     ps.x = bsxfun(@(x,y) x + cos(y),ps.x,ps.phy_heading);
        %     ps.y = bsxfun(@(x,y) x + sin(y),ps.y,ps.phy_heading);
        %     ps.x = ps.x + cos(ps.mag_heading-pi/2);     % when heading shifted input date
        %     ps.y = ps.y + sin(ps.mag_heading-pi/2);
            hN = .03;      % heading noise value
            if i>1
        %         Halpha = pi-est(i-1,3);
                Halpha = 0;
                ps.mag_heading = ps.mag_heading+random('normal',Halpha,hN,n,1);     % candidate, .08
            else
        %     ps.mag_heading = ps.mag_heading+euler(i,3);
                ps.mag_heading = ps.mag_heading+random('normal',0,hN,n,1);     % candidate, .08
            end

            ps.sl = .7 + random('normal',0,.5,n,1);

            mu = [0,0];
        %     sigma = [0.10409786, 0.13461109; 0.13461109, 0.29744705];
        %     sigma = [1.43724175 -0.93884837;-0.93884837  0.74606608];
            sigma = [0.02765426 0;0  0.04181993];
            mvnRand = mvnrnd(mu,sigma,n);
        %     ps.x = ps.x + cos(ps.mag_heading+euler(i,3)+yaw_offset).*ps.sl+ mvnRand(:,1);
        %     ps.y = ps.y + sin(ps.mag_heading+euler(i,3)+yaw_offset).*ps.sl+ mvnRand(:,2);

            ps.x = ps.x + cos(ps.mag_heading+euler(i,3)).*ps.sl+ random('Uniform',-.1,.1,n,1);
            ps.y = ps.y + sin(ps.mag_heading+euler(i,3)).*ps.sl+ random('Uniform',-.1,.1,n,1);

        %     ps.x = ps.x + cos(ps.mag_heading+euler(i,3))*sl;
        %     ps.y = ps.y + sin(ps.mag_heading+euler(i,3))*sl;
        %     ps.x = ps.x + cos(ps.mag_heading+euler(i,3))*sl + random('Uniform',-1,1,n,1);
        %     ps.y = ps.y + sin(ps.mag_heading+euler(i,3))*sl + random('Uniform',-1,1,n,1);
        %     ps.x = ps.x + cos(ps.mag_heading+d_heading(i))*sl + random('Uniform',-1,1,n,1);
        %     ps.y = ps.y + sin(ps.mag_heading+d_heading(i))*sl + random('Uniform',-1,1,n,1);
        %     ps.x = ps.x + ps.stlng.*cos(ps.heading);
        %     ps.y = ps.y + ps.stlng.*sin(ps.heading);

            % ================ UPDATE
            % 1. find (geo-locational) nearest learning data
            [phy_dist,I] = pdist2([lm.x,lm.y],[ps.x,ps.y],'euclidean','Smallest',1);
            % 2. calculate Rotated magnetic field data and magnetic distance

        %     R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]/(rotMat(:,:,i))),0,...
        %         'UniformOutput',false);
        %     rotatedMag = cell2mat(cellfun(@(x)((x*tM(i,:)')'),R,'UniformOutput',false));

        %     (1) UPDATE FUNCTION candidate 
        %     R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]*(rotMat(:,:,i))),ps.mag_heading,...
        %         'UniformOutput',false); 
        %     rotatedMag = cell2mat(cellfun(@(x)(x.'*tM(i,:)')',R,'UniformOutput',false));
        %     (2) UPDATE FUNCTION candidate  
        %     R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]/(rotMat(:,:,i))),ps.mag_heading,...
        %         'UniformOutput',false); 
        %     rotatedMag = cell2mat(cellfun(@(x)(x*tM(i,:)')',R,'UniformOutput',false));
        %     (3) UPDATE FUNCTION candidate  
        %     R = arrayfun(@(x)(rotMat(:,:,i)*[cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]),...
        %         -ps.mag_heading,'UniformOutput',false); 
        %     rotatedMag = cell2mat(cellfun(@(x)(x.'*tM(i,:)')',R,'UniformOutput',false));

        %     R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1])*rotMat(:,:,i),...
        %         ps.mag_heading,'UniformOutput',false); 
        %     rotatedMag = cell2mat(cellfun(@(x)(x.'*tM(i,:)')',R,'UniformOutput',false));

        %     R = arrayfun(@(x)([cos(x) -sin(x) 0;sin(x) cos(x) 0;0 0 1]*(rotMat(:,:,i).')),ps.mag_heading,...
        %         'UniformOutput',false); 
        %     rotatedMag = cell2mat(cellfun(@(x)(x*tM(i,:)')',R,'UniformOutput',false));
            % TODO: ???validated?
        %     rotatedMag2 = getHeadingRotatedVector(ps.mag_heading, tM(i,:), rotMat(:,:,i));
        %     if ~isequal(rotatedMag,rotatedMag2)
        %         disp('!!')
        %     end
        %     rotatedMag = rotatedMag2;

        %     R = arrayfun(@(x) eulerAnglesToRotation3d(x,-euler(1,2),-euler(1,1))...
        %         ,-euler(1,3)-ps.mag_heading,'UniformOutput',false);
        %     rotatedMag = cell2mat(cellfun(@(x)(x*tM(i,:)')',R,'UniformOutput',false));    
        %     R = arrayfun(@(x) euler2rotMat(-euler(1,1),-euler(1,2),x),...
        %         -euler(1,3)-ps.mag_heading,'UniformOutput',false);
        %     rotatedMag = cell2mat(cellfun(@(x)(x*tM(i,:)')',R,'UniformOutput',false));
        %     R = arrayfun(@(x) euler2rotMat(euler(1,1),euler(1,2),x),...
        %         euler(1,3)+ps.mag_heading,'UniformOutput',false);
        %     rotatedMag = cell2mat(cellfun(@(x)(x.'*tM(i,:)')',R,'UniformOutput',false));

            % EUCLIDEAN
            % TODO: may be more optimizable (DONE?maybe)
        %     mag_dist = diag(pdist2(rotatedMag,lM(I,:),'euclidean'));    
        %     mag_dist = bsxfun(@(x,y) pdist([x;y]), rotatedMag,lM(I,:));
        %     mag_dist = sqrt(sum((rotatedMag-lM(I,:)).^2,2));

            % ILoA, line 256, 260, 276
            rotatedMag2 = getHeadingRotatedVector(ps.mag_heading, tM(i,:), rotMat(:,:,i));
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
            ps.prob = 1./(mag_dist);
            in = isinterior(shp,ps.x,ps.y);
        %     in = inShape(shp,[ps.x,ps.y]);
            ps.prob(~in) = 0;
            ps.prob(phy_dist'>3) = 0;

            if sum(ps.prob) == 0
                rand_idx = randi(length(lm.x),n,1);
                ps.x = lm.x(rand_idx);
                ps.y = lm.y(rand_idx);
                ps.mag_heading = random('Uniform', 0,2*pi,n,1);
        %         ps.phy_heading = random('Uniform', 0,2*pi,n,1);
                ps.prob = ones(n,1)*(1/n);
            else
                ps.prob = ps.prob./sum(ps.prob);
            end

            % ================ RESAMPLE
            resample_idx = randsample(1:n,n,true,ps.prob);
        %     phy_move_noise_range = 2;
        %     ps.x = ps.x(resample_idx) + phy_move_noise_range*rand(n,1) - phy_move_noise_range/2;
        %     ps.y = ps.y(resample_idx) + phy_move_noise_range*rand(n,1) - phy_move_noise_range/2;
        %     ps.x = ps.x(resample_idx) + random('normal',0,.5,n,1);
        %     ps.y = ps.y(resample_idx) + random('normal',0,.5,n,1);
            ps.x = ps.x(resample_idx);
            ps.y = ps.y(resample_idx);
            ps.mag_heading = ps.mag_heading(resample_idx);
            ps.sl = ps.sl(resample_idx);

        %     if std_euler(i) > 0.03
        %         ps.mag_heading = ps.mag_heading(resample_idx)+random('normal',0,std_euler(i),n,1);
        %     else
        %         ps.mag_heading = ps.mag_heading(resample_idx);
        %     end
        %     ps.mag_heading = ps.mag_heading(resample_idx)+random('normal',0,heading_noise,n,1);    
        %     ps.mag_heading = ps.mag_heading(resample_idx);
        %     ps.phy_heading = ps.phy_heading(resample_idx)+random('normal',0,.001,n,1);

        %     ps.heading = ps.heading(resample_idx)+random('Uniform', -pi/10,pi/10,n,1);
        %     ps.stlng = ps.stlng(resample_idx) + random('normal',0,1,n,1);
        %     ps.prob = init_prob;

            set(h_ps,'XData',ps.x,'YData',ps.y)                             % ps result
            set(h_pm,'XData',mean(ps.x),'YData',mean(ps.y));
        %     set(h_gt,'XData',data2.coord_x(i),'YData',data2.coord_y(i))     % ground truth
            drawnow

        %     est(i,:) = [mean(ps.x),mean(ps.y),mean(angdiff(ps.mag_heading,0))];
        %     est(i,:) = [mean(ps.x),mean(ps.y),mean(abs(angdiff(ps.mag_heading,0)))];        %x,y,heading
        %     est(i,:) = [mean(ps.x),mean(ps.y),circ_mean(ps.mag_heading-pi)]; 
            est(i,:) = [mean(ps.x),mean(ps.y),circ_mean(ps.mag_heading)]; 
            eta(i,:) = wrapTo2Pi(ps.mag_heading'); 

        %     err(i,:) = pdist2([data2.coord_x(i),data2.coord_y(i)],[ps.x,ps.y]);
            AD(i,:) = pdist2([mean(ps.x),mean(ps.y)],[ps.x,ps.y]);
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
hold on
for i=1:4
    ecdf(vertcat(errs{:,i})/pi*180)
end
legend('85','45','0','20')