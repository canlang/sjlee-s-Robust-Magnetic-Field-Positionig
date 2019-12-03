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

map = magmap_construction('mats',.5);
plot(map(:,1),map(:,2),'x')