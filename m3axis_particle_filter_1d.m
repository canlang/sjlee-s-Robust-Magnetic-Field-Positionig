clear; close all

video_flag = false;
with_histogram = false;

mag2 = readtable('db_20171020_15_48_22/magnetic.csv');
mag4 = readtable('db_20171020_15_51_04/magnetic.csv');

% rotZrad = deg2rad(180);
% rotZ = [cos(rotZrad), -sin(rotZrad), 0;
%     sin(rotZrad), cos(rotZrad), 0;
%     0, 0, 1];
lM = [smooth(mag2.x,10), smooth(mag2.y,10), smooth(mag2.z,10)];
tM = [smooth(mag4.x,10),smooth(mag4.y,10),-smooth(mag4.z,10)];

n = 500;
x = random('Uniform', 1,length(lM),n,1);
%%% DATA INFO. ps: [location, probability, heading(forward & backward)]
ps = [x,ones(n,1)*(1/n), randi([0 1], n,1)];


% Video save
if video_flag
    v = VideoWriter('m3axis_1d_pf.mp4');
    v.FrameRate = 10;
    v.Quality = 100;
    open(v);
end

% co = [215,25,28;
%     253,174,97;
%     44,123,182]/255;

co = [55,126,184;
77,175,74;
255,127,0]/255;
set(0,'DefaultAxesColorOrder', co,...
      'DefaultAxesLineStyleOrder','-|:|--')

switch with_histogram
    case true
        subplot(211)
end

plot(lM)
legend()
% plot(lM(:,1),'LineWidth',2)     % learning data plot
% plot(lM(:,2:3))
xlabel('Data index: sampling freq. = 1.87 (ms)')
ylabel('Magnetic field (\muT)')
ylim([-70 70])
xlim([0 3000])

% set(groot,'defaultAxesColorOrder','remove')
% co = [253,174,97;
% 255,255,191;
% 44,123,182;]/255;
% set(0,'DefaultAxesColorOrder', co)
% ax = axes('ColorOrder',co,'NextPlot','replacechildren');
hold on
plot(-tM)   % test data plot
h_ps = plot(ps(:,1), ones(n,1)*60, '.','MarkerSize', 12,'color',[215,25,28]/255);    % particle plot
% h_ps = plot(ps(:,1), zeros(length(ps)), 'r.', 'MarkerSize', 12);    % particle plot
h_l = refline([0 -tM(1)]);  % ref line plot
h_l.Color = [171,217,233]/255;
h_input_idx = line([0,0],ylim,'Color',[171,217,233]/255,'LineStyle','--');
% legend({'x comp. as Learning','y comp. as Learning','z comp. as Learning',...
%     'x comp. as Testing','y comp. as Testing','z comp. as Testing',...
%     'Particle (prob.)','x comp. input reference'},'Location','southeast','NumColumns',3)
legend({'|m_\alpha| (learning)','|m_\beta| (learning)','|m_\gamma| (learning)',...
    '|b_\alpha| (testing)','|b_\beta| (testing)','|b_\gamma| (testing)',...
    '\bf{\omega} of particles','Current reference'},'Location','southeast',...
    'FontName','Source Code Pro','NumColumns',3)
% grid on
hold off    

switch with_histogram
    case true        
        subplot(212)
        h_hist = histogram(ps(:,2),0:0.01:1);
        ylabel('Bin count')
        xlabel('Probability')
end

set(gcf,'units','points','position',[800,500,800,400])
sdf(gcf,'sj')
%%
ST = 5;
% for i = 1:200
for i = 1:ST:length(tM)
    % predict
    CON = ps(:,3)==0;
    
    ps(CON,1) = ps(CON,1)-ST;
    ps(~CON,1) = ps(~CON,1)+ST;
    
    
    location = round(ps(:,1));
    % UPDATE v1
    for j = 1:length(ps)
        if location(j)>0 && location(j)<length(lM)
            if ps(j,3) == 1
                rot = -1;
            else
                rot = 1;
            end
            D = pdist([rot*tM(i,:);lM(location(j),:)]);
            prob = 1/D;
%             prob = 1/abs(rot*tM(j) - lM(location(j)));
            if prob ~= inf
%                 ps(j,2) = prob*ps(j,2);
                ps(j,2) = prob;
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
    
    % UPDATE v2
%     rotatedlM = lM(location,:);
%     ps(:,2) = 1./pdist2(tM(i,:),rotatedlM);
%     ps(:,2) = ps(:,2)/sum(ps(:,2)); % normalize PROP.
    
    % RESAMPLING v1
%     RCON = ps(:,2) <= 0.001;
%     nResample = nnz(RCON);
%     restProb = sum(ps(RCON,2));
%     ps(RCON,:) = [];
%     if nResample > 0
% %         newLoc = random('Uniform', 1,length(mag2.time),nResample,1);
% %         newPs = [newLoc, restProb*ones(nResample,1), randi([0 1], nResample,1)];
%         newLoc = randsample(ps(:,1),nResample, true, ps(:,2))+10*rand(nResample,1);
%         newPs = [newLoc, 0.0001*ones(nResample,1), randi([0 1], nResample,1)];
%         ps = [ps;newPs];
%     end
    
    % RESAMPLING v2
    resample_idx = randsample(1:n,n,true,ps(:,2));
    ps(:,1) = ps(resample_idx,1) + random('normal',0,1,n,1);
    ps(:,3) = ps(resample_idx,3);
    
    set(h_ps,'XData', ps(:,1))
    set(h_l,'YData', [-tM(i) -tM(i)])
    h_input_idx.XData = [i i];
%     subplot(212)
    switch with_histogram
        case true       
            set(h_hist,'Data', ps(:,2))
    end
    drawnow limitrate
    
    if video_flag
        frame = getframe(gcf);
        writeVideo(v,frame);
    end

    if i >= 15
        break
    end
end
if video_flag close(v);end

tightfig
print -depsc2 eps/3axis_1d_ambiguity15.eps