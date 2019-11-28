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

mu = [0 0];
% sigma = [1 1.5; 1.5 3];
% sigma = [1 1; 1 10];
sigma = [0.1,0;0,0.01];
rng('default')  % For reproducibility
R = mvnrnd(mu,sigma,10000);
plot(R(:,1),R(:,2),'+')
axis equal