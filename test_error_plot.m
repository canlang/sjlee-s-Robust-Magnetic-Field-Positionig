err_std = std(errMat,0,2);
err_mean = mean(errMat,2);
errorbar(nParticleCandidate,err_mean,err_std)