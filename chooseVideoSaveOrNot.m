function video_flag= chooseVideoSaveOrNot
answer = questdlg('Would you like save the video result?', ...
	'Record the result', ...
	'Yes','No','No');
switch answer
    case 'Yes'
        disp([answer ', video file would be create or overwrite.'])
        video_flag = 1;
    case 'No'
        disp([answer ', just see the result.'])
        video_flag = 0;
end