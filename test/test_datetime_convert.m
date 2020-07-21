%%
% datetime(res_data.Time(step_idx),'InputFormat', 'yyyy-MM-dd HH:mm:ss.SSS +0900')
% posixtime(res_data.Time(step_idx))
% datetime(res_data.Time,'ConvertFrom','posixtime','TimeZone','local')
% datetime(res_data.Time(step_idx),'ConvertFrom','posixtime','TimeZone','')

% datetime(res_data.Time(step_idx),'ConvertFrom','posixtime','TimeZone','America/Los_Angeles')
% datetime(res_data.Time(step_idx),'ConvertFrom','epochtime','Epoch','2020-01-01')
% time_sec = seconds(res_data.Time);
% time_sec.Format = 'hh:mm'; 

% a = datetime(rawdata.mag(:,1)*1e-3,'ConvertFrom','posixtime','TimeZone','Asia/Seoul')
% b = datetime(rawdata.mag(:,2)*1e-6,'ConvertFrom','posixtime','TimeZone','Asia/Seoul')

% https://developer.android.com/reference/android/os/SystemClock
% second row of rawdata is elapse since system start, therefore, there is
% no meaning. BUT, this is because posix time (first row) have same (not
% unique) value

%%
% d = datetime('now');
% posixtime(d)