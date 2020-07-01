%% local functions
function data = load_rawdata(datapath,varargin)
if nargin < 1
    print('hi')
end
if nargin == 1
    acc = csvread(fullfile(datapath,'Accelerometer.csv'),1,0);
    gyr = csvread(fullfile(datapath,'Gyroscope.csv'),1,0);
    mag = csvread(fullfile(datapath,'MagneticField.csv'),1,0);
%     acc(:,2) = acc(:,2)*10^3;
%     gyr(:,2) = gyr(:,2)*10^3;
%     mag(:,2) = mag(:,2)*10^3;
    data.acc = acc;
    data.gyr = gyr;
    data.mag = mag;
%     data.acc_time = datetime(acc(:,1)/10^3,'convertfrom','posixtime','TimeZone','Asia/Seoul');
%     data.gyr_time = datetime(gyr(:,1)/10^3,'convertfrom','posixtime','TimeZone','Asia/Seoul');
%     data.mag_time = datetime(mag(:,1)/10^3,'convertfrom','posixtime','TimeZone','Asia/Seoul');

%     data.acc_norm = vecnorm(acc(:,3:5),2,2);
%     data.acc_norm = data.acc_norm - mean(data.acc_norm);
end

if (nargin == 2 && strcmp(varargin{1},'iPhone'))
%     filename = fullfile(datapath.folder,datapath.name);
    filename = datapath;
    
%     M = readmatrix(filename,'Delimiter',{','});
%     M = unique(M,'rows');
%     acc = M(:,40:42)+M(:,48:50);     %TODO : second '2' column is not meaning, so have to re-designed
%     gyr = M(:,37:39);
%     mag = M(:,51:53);
%     times = [M(:,14)*10^3,M(:,14)*10^6];    
%     data.acc = [times,acc];
%     data.gyr = [times,gyr];
%     data.mag = [times,mag];
% 
%     times = M(1,14);    
%     tacc = M(:,33)-M(1,33)+times;
%     tgyr = M(:,33)-M(1,33)+times;
%     tmag = M(:,33)-M(1,33)+times;
% %     data.acc = [tacc*10^3,tacc*10^6,acc*9.8];
% %     data.gyr = [tgyr*10^3,tgyr*10^6,gyr];
% %     data.mag = [tmag*10^3,tmag*10^6,mag];
%     data.acc = unique([tacc*10^3,tacc*10^6,acc*9.8],'rows');
%     data.gyr = unique([tgyr*10^3,tgyr*10^6,gyr],'rows');
%     data.mag = unique([tmag*10^3,tmag*10^6,mag],'rows');
     
    T = readtable(filename,'Delimiter',',','PreserveVariableNames',true);
    dt = datetime(T{:,1},'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS +0900');
    pt = posixtime(dt);
    
    acc = [T.("accelerometerAccelerationX(G)"),T.("accelerometerAccelerationY(G)"),T.("accelerometerAccelerationZ(G)")];
%     acc = [T.("accelerometerAccelerationX(G)"),T.("accelerometerAccelerationY(G)"),T.("accelerometerAccelerationZ(G)")]...
%         +[T.("motionGravityX(G)"),T.("motionGravityY(G)"),T.("motionGravityZ(G)")];
    gyr = [T.("gyroRotationX(rad/s)"),T.("gyroRotationY(rad/s)"),T.("gyroRotationZ(rad/s)")];
%     mag = [T.("motionMagneticFieldX(µT)"),T.("motionMagneticFieldY(µT)"),T.("motionMagneticFieldZ(µT)")];
    mag = [T.("locationHeadingX(µT)"),T.("locationHeadingY(µT)"),T.("locationHeadingZ(µT)")];
    
%     acc = [T.("motionUserAccelerationX(G)"),T.("motionUserAccelerationY(G)"),T.("motionUserAccelerationZ(G)")];
%     gyr = [T.("motionRotationRateX(rad/s)"),T.("motionRotationRateY(rad/s)"),T.("motionRotationRateZ(rad/s)")];
%     mag = [T.("motionMagneticFieldX(µT)"),T.("motionMagneticFieldY(µT)"),T.("motionMagneticFieldZ(µT)")];
    
%     acc = T{:,40:42}+T{:,48:50};     %TODO : second '2' column is not meaning, so have to re-designed
%     26,37
%     gyr = -T{:,26:28};
%     gyr = -T{:,37:39};
%     30
%     mag = T{:,51:53};    
    
%     data.acc = [pt*10^3,pt,acc*9.8];
%     data.gyr = [pt*10^3,pt,gyr];
%     data.mag = [pt*10^3,pt,mag];
%     data.acc = [pt*10^3,pt*10^6,acc*9.8];
%     data.gyr = [pt*10^3,pt*10^6,gyr];
%     data.mag = [pt*10^3,pt*10^6,mag];
    data.acc = [pt*10^3,pt*10^9,acc*9.8];
    data.gyr = [pt*10^3,pt*10^9,gyr];
    data.mag = [pt*10^3,pt*10^9,mag];
    
    data.quaternion = [T.("motionQuaternionW(R)"),T.("motionQuaternionX(R)"),T.("motionQuaternionY(R)")...
        ,T.("motionQuaternionZ(R)")];
    data.euler = [T.("motionRoll(rad)"),T.("motionPitch(rad)"),T.("motionYaw(rad)")];
%     ];
end

data.acc_time = datetime(data.acc(:,1)/10^3,'convertfrom','posixtime','TimeZone','Asia/Seoul');
data.gyr_time = datetime(data.gyr(:,1)/10^3,'convertfrom','posixtime','TimeZone','Asia/Seoul');
data.mag_time = datetime(data.mag(:,1)/10^3,'convertfrom','posixtime','TimeZone','Asia/Seoul');

data.acc_norm = vecnorm(data.acc(:,3:5),2,2);
data.acc_norm = data.acc_norm - mean(data.acc_norm);
end