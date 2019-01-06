clear; close all
mag2 = readtable('db_20171020_15_48_22/magnetic.csv');
mag4 = readtable('db_20171020_15_51_04/magnetic.csv');

% rotZrad = deg2rad(180);
% rotZ = [cos(rotZrad), -sin(rotZrad), 0;
%     sin(rotZrad), cos(rotZrad), 0;
%     0, 0, 1];
lM = [smooth(mag2.x,10), smooth(mag2.y,10), smooth(mag2.z,10)];
tM = [smooth(mag4.x,10),smooth(mag4.y,10),-smooth(mag4.z,10)];

n = 400;
x = random('Uniform', 1,length(lM),n,1);
ps = [x,ones(n,1)*(1/n), randi([0 1], n,1)];
% (DATA INFO) ps: [location, probability, heading]

% Video save
% v = VideoWriter('m3axis_1d_pf.mp4');
% v.FrameRate = 10;
% v.Quality = 100;
% open(v);

% learning data plot
subplot(211)
plot(lM(:,1),'LineWidth',2)
hold on
plot(lM(:,2:3))
set(gcf,'units','points','position',[800,500,800,500])
ylim([-50 80])
xlim([0 3000])

% test data plot
plot(-tM)
% particle plot
h_ps = plot(ps(:,1), zeros(length(ps)), 'r.', 'MarkerSize', 12);
% ref line plot
h_l = refline([0 -tM(1)]);
% h_l.Color = 'r';
h_l.Color = [.5 0 .9];
hold off    
subplot(212)
hist(ps(:,2))
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
            if ps(j,3) == 1
                rot = -1;
            else
                rot = 1;
            end
            D = pdist([rot*tM(j,:);lM(location(j),:)]);
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
    
    RCON = ps(:,2) <= 0.001;
    nResample = nnz(RCON);
    restProb = sum(ps(RCON,2));
    ps(RCON,:) = [];
    if nResample > 0
%         newLoc = random('Uniform', 1,length(mag2.time),nResample,1);
        newLoc = randsample(ps(:,1),nResample, true, ps(:,2))+10*rand(nResample,1);
%         newPs = [newLoc, restProb*ones(nResample,1), randi([0 1], nResample,1)];
        newPs = [newLoc, 0.0001*ones(nResample,1), randi([0 1], nResample,1)];
        ps = [ps;newPs];
    end
    
    
    set(h_ps,'XData', ps(:,1))
    set(h_l,'YData', [-tM(i) -tM(i)])
    subplot(212)
    hist(ps(:,2))
    drawnow limitrate
    
    % Video save
%     frame = getframe(gcf);
%     writeVideo(v,frame);
end
close(v)