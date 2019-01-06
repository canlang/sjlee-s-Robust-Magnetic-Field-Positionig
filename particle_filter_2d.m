clear; close all
mag2 = readtable('db_20171020_15_48_22/magnetic.csv');
mag4 = readtable('db_20171020_15_51_04/magnetic.csv');

% rotZrad = deg2rad(180);
% rotZ = [cos(rotZrad), -sin(rotZrad), 0;
%     sin(rotZrad), cos(rotZrad), 0;
%     0, 0, 1];
lM = smooth(mag2.x,10);
tM = smooth(mag4.x,10);

n = 300;
x = random('Uniform', 1,length(mag2.time),n,1);
ps = [x,ones(n,1)*(1/n), randi([0 1], n,1)];
% ps: [location, probability, heading]

% for i = 1:200
for i = 1:5:length(tM)
    % predict
    CON = ps(:,3)==0;
    ps(CON,1) = ps(CON,1)-1;
    ps(~CON,1) = ps(~CON,1)+1;
    
    % update
    location = round(ps(:,1));    
    for j = 1:length(ps)
        if location(j)>0 && location(j)<length(lM)
            if ps(j,3) == 1
                rot = -1;
            else
                rot = 1;
            end
            prob = 1/abs(rot*tM(j) - lM(location(j)));
            if prob ~= inf
                ps(j,2) = prob;
            else
                % prob
                ps(j,2) = 50;
            end
            
%             pd = makedist('Normal', lM(location(j)), 1);
%             ps(j,2) = pdf(pd, -tM(j));
        else
            ps(j,2) = 0;
        end
    end
    ps(:,2) = ps(:,2)/sum(ps(:,2));
    
    nResample = nnz(ps(:,2) < 0.005);
    ps(ps(:,2)<0.005,:) = [];
    if nResample > 0
%         newLoc = random('Uniform', 1,length(mag2.time),nResample,1);
        newLoc = randsample(ps(:,1),nResample, true, ps(:,2))+10*rand(nResample,1);
        newPs = [newLoc, zeros(nResample,1), randi([0 1], nResample,1)];
        ps = [ps;newPs];
    end
    
    
    plot(lM)
    hold on
    plot(-tM)
    ylim([-50 80])
    
    plot(ps(:,1), zeros(length(ps)), 'r.')
%     plot(i,-tM(i),'o')
    hold off
    
%     delete(h)
    h = refline([0 -tM(i)]);
    h.Color = 'r';
%     h = refline([inf -mag4.x(i)]);
%     h.Color = 'r';
    drawnow limitrate
    
end