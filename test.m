% [inputs,targets] = cancer_dataset;
% hiddenLayerSize = 10;
% net = patternnet(hiddenLayerSize);
% 
% net.divideParam.trainRatio = 70/100;
% net.divideParam.valRatio   = 15/100;
% net.divideParam.testRatio  = 15/100;
% 
% 
% [net,tr] = train(net,inputs,targets);
% 
% outputs = net(inputs);
% errors = gsubtract(targets,outputs);
% performance = perform(net,targets,outputs)
% 
% tInd = tr.testInd;
% tstOutputs = net(inputs(:,tInd));
% tstPerform = perform(net,targets(:,tInd),tstOutputs)

% ---------------------------------------------------------------------------
% map = magmap_construction('mats',.5);
% plot(map(:,1),map(:,2),'x')
% ---------------------------------------------------------------------------

err_cnt = 0;
for i = 1:1000
    try 
        a = rand(3); 
        c = cov(a) + .0001 * eye(3); 
        m = mean(a); 
        mvnpdf(a, m, c);
    catch me
        err_cnt = err_cnt + 1;
    end
end
%%
close all
n = 1000;
mu = [0,0];
% sigma = [0.10409786, 0.13461109; 0.13461109, 0.29744705];
% sigma = [0.27791309 0.22856106;0.22856106 0.16106075];
sigma = [2.97925885  0.02243997;0.02243997  0.98844656];
eig(sigma)
sigma - sigma.'
S = sigma' * sigma
S2 = (sigma + sigma')/2
%%
mvnRand = mvnrnd(mu,sigma,n);
plot(mvnRand(:,1),mvnRand(:,2),'+')
axis equal
sdf(gcf,'sj')