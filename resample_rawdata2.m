function data = resample_rawdata2(rawdata,rate)
% Var1 : 3 vectors of Accelerometer, Var2: norm vector of Acc.

T_acc = timetable(seconds(rawdata.acc(:,1)/1e3),rawdata.acc(:,3:5),rawdata.acc_norm,...
    'VariableNames',{'acc','acc_norm'});
% Var1 : 3 vectors of gyroscope
T_gyr = timetable(seconds(rawdata.gyr(:,1)/1e3),rawdata.gyr(:,3:5),...
    'VariableNames',{'gyr'});
T_mag = timetable(seconds(rawdata.mag(:,1)/1e3),rawdata.mag(:,3:5),...
    'VariableNames',{'mag'});

% TODO: retime?
% https://kr.mathworks.com/help/matlab/matlab_prog/clean-timetable-with-missing-duplicate-or-irregular-times.html

% T_acc = sortrows(T_acc);
% T_gyr = sortrows(T_gyr);
% T_mag = sortrows(T_mag);
[~,ia] = unique(T_acc.Time);
T_acc = T_acc(ia,:);
[~,ia] = unique(T_gyr.Time);
T_gyr = T_gyr(ia,:);
[~,ia] = unique(T_mag.Time);
T_mag = T_mag(ia,:);

TT = synchronize(T_acc,T_gyr,T_mag,'regular','linear','TimeStep',seconds(rate));

data.Accelerometer = TT.acc;
data.Gyroscope = TT.gyr;
% data.Gyroscope = TT.gyr*(180/pi);
data.Magnetometer = TT.mag;
data.acc_norm = TT.acc_norm;

data.Time = seconds(TT.Time(:));
% data.Time = seconds(TT.Time(:)-(TT.Time(1)));
data.Rate = median(diff(data.Time)); % cal sample rate
if isequaln(data.Rate, rate)
    fprintf(1,'Trying to set %.2f rate, but %.2f setted.\n',data.Rate,rate);
end
end