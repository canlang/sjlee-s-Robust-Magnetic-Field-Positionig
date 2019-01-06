% 1-axis mean magnitude value of geo-magnetic field vector
clear; close all
mag2 = readtable('db_20171020_15_48_22/magnetic.csv');
mag4 = readtable('db_20171020_15_51_04/magnetic.csv');

lM = [smooth(mag2.x,10), smooth(mag2.y,10), smooth(mag2.z,10)];
tM = [smooth(mag4.x,10),smooth(mag4.y,10),-smooth(mag4.z,10)];


% 1-1. input as x component of magnetic field
lM = lM(:,1);
tM = abs(tM(:,1));

% % 1-2. calculate magnitude value (L2norm)
% lM = vecnorm(lM')';
% tM = vecnorm(tM')';

n = 300;
x = random('Uniform', 1,length(lM),n,1);
ps = [x,ones(n,1)*(1/n), randi([0 1], n,1)];
% (DATA INFO) ps: [location, probability, heading]


% v = VideoWriter('m1axis_1d_pf_sis.mp4','MPEG-4');
% v.FrameRate = 10;
% v.Quality = 100;
% open(v);

% learning data plot
% subplot(211)
plot(lM,'LineWidth',2)
set(gcf,'units','points','position',[800,500,800,500])
ylim([-50 80])
xlim([0 3000])
hold on
% test data plot
plot(tM)
% particle plot
h_ps = plot(ps(:,1), zeros(length(ps)), 'r.','MarkerSize',12);
% ref line plot
h_l = refline([0 tM(1)]);
h_l.Color = 'r';
hold off    
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
            D = pdist([tM(i);lM(location(j))]);
            prob = 1/D;
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
    ps(:,2) = ps(:,2)/sum(ps(:,2));
    if sum(ps(:,2)) ~= 1
        ps(:,2) = 1/length(ps);
    end
    
    noise_range = 3;
    resampleIdx = randsample(1:n,n, true, ps(:,2));
    ps(:,[1 3]) = ps(resampleIdx,[1 3]);
    ps(:,1) = ps(:,1) + noise_range*rand(n,1) - noise_range/2;
    ps(:,2) = ones(length(ps),1)/n;
    
    
    set(h_ps,'XData', ps(:,1),'YData',normrnd(0,2,[n,1]))
    set(h_l,'YData', [tM(i) tM(i)])
%     subplot(212)
%     hist(ps(:,2))
    drawnow limitrate
%     frame = getframe(gcf);
%     writeVideo(v,frame);
%     break
    if i >= 10
        break
    end
end
% close(v);

xlabel('Data index')
ylabel("x component (\muT)")
sdf(gcf,'sj2')
tightfig(gcf)