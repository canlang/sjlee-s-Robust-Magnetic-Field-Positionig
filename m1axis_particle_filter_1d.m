% 1-axis mean magnitude value of geo-magnetic field vector
clear; close all

video_flag = false;
mag2 = readtable('db_20171020_15_48_22/magnetic.csv');
mag4 = readtable('db_20171020_15_51_04/magnetic.csv');
% datetime(mag2.time/1000,'ConvertFrom','posixtime','TimeZone','Asia/Seoul')
lM = [smooth(mag2.x,10), smooth(mag2.y,10), smooth(mag2.z,10)];
tM = [smooth(mag4.x,10),smooth(mag4.y,10),-smooth(mag4.z,10)];

% 1-1. input as x component of magnetic field
% lM = lM(:,1);
% tM = abs(tM(:,1));

% % 1-2. calculate magnitude value (L2norm)
lM = vecnorm(lM')';
tM = vecnorm(tM')';

n = 500;
x = random('Uniform', 1,length(lM),n,1);
ps = [x,ones(n,1)*(1/n), randi([0 1], n,1)];
% (DATA INFO) ps: [location, probability, heading]

switch video_flag
    case true
        v = VideoWriter('m1axis_1d_pf.mp4');
        v.FrameRate = 10;
        v.Quality = 100;
        open(v);
end
% learning data plot
% subplot(211)
plot(lM,'color',[44,123,182]/255)
xlabel('Data index: sampling freq. = 1.87 (ms)')
ylabel('Magnetic field (\muT)')

ylim([0 90])

% test data plot
hold on
plot(tM,':','color',[44,123,182]/255)
h_ps = plot(ps(:,1), ones(n,1)*80, '.','MarkerSize',12,'color',[215,25,28]/255);      % particle plot
h_l = refline([0 tM(1)]);       % ref line plot
h_l.Color = [171,217,233]/255;
% h_input_idx = line([0,0],ylim,'Color',[.5 0 .9],'LineStyle','--');
hold off 

legend({'|\bf{m}| (learning)','|\bf{b}| (testing)',...
    '\bf{\omega} of particles','Current reference'},'Location','southeast',...
    'FontName','Source Code Pro')
set(gcf,'units','points','position',[800,400,800,400])
sdf(gcf,'sj')

% subplot(212)
% hist(ps(:,2))
%%
ST = 5;
% for i = 1:200
for i = 1:ST:length(tM)
    % predict
    CON = ps(:,3)==0;
    
    ps(CON,1) = ps(CON,1)-ST;
    ps(~CON,1) = ps(~CON,1)+ST;
    
    % update
    location = round(ps(:,1));    
    for j = 1:length(ps)
        if location(j)>0 && location(j)<length(lM)
%             if ps(j,3) == 1
%                 rot = -1;
%             else
%                 rot = 1;
%             end
            D = pdist([tM(i);lM(location(j))]);
            prob = 1/D;
%             prob = 1/abs(rot*tM(j) - lM(location(j)));
            if prob ~= inf
                ps(j,2) = prob*ps(j,2);
%                 ps(j,2) = prob;
            else
%                 prob
                ps(j,2) = 1;
            end
            
%             pd = makedist('Normal', lM(location(j)), 1);
%             ps(j,2) = pdf(pd, -tM(j));
        else
            ps(j,2) = 0;
        end
    end
    ps(:,2) = ps(:,2)/sum(ps(:,2));
    
    % RESAMPLING v1
%     RCON = ps(:,2) <= 0.001;
%     nResample = nnz(RCON);
%     restProb = sum(ps(RCON,2));
%     ps(RCON,:) = [];
%     if nResample > 0
% %         newLoc = random('Uniform', 1,length(mag2.time),nResample,1);
%         newLoc = randsample(ps(:,1),nResample, true, ps(:,2))+50*rand(nResample,1);
%         newPs = [newLoc, restProb*ones(nResample,1), randi([0 1], nResample,1)];        
%         ps = [ps;newPs];
%     end
    
    % RESAMPLING v2
    resample_idx = randsample(1:n,n,true,ps(:,2));
    ps(:,1) = ps(resample_idx,1) + random('normal',0,2,n,1);
    ps(:,3) = ps(resample_idx,3);
    
    
    set(h_ps,'XData', ps(:,1))
    set(h_l,'YData', [tM(i) tM(i)])
    h_input_idx.XData = [i i];
    drawnow limitrate
    
    switch video_flag
        case true
            frame = getframe(gcf);
            writeVideo(v,frame);
    end
    
    if i >= 13
        break
    end
end
switch video_flag
    case true
        close(v);
end

tightfig
print -depsc2 eps/1axis_1d_ambiguity15.eps

% set(0,'DefaultAxesColorOrder',brewermap(2,'Set2'))
% tightfig(gcf)
% set(0,'DefaultAxesColorOrder','remove')