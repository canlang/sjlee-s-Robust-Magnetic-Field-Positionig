
Nloop = 10;
Nfailure = 0;
Nsuccess = 0;
errs = cell(Nloop,1);

intp = [.1, .2, .3, .5, .8, 1.0, 1.2];
for j=1:length(intp)
    for i=1:Nloop
        close all
        err = ILoA('KI-1F','MATE20pro',3,intp(j),false);
        %%    
        if sum(err<3)>length(err)/2
            % disp converged!
            Nsuccess = Nsuccess+1;
            errs{i} = err;
        else
            Nfailure = Nfailure+1;
            % disp fail!
        end
    end
    all_errs = errs{:};
    MED = mean(rmoutliers(all_errs));
    fprintf('Convergence rate: %.1f, MED: %.2f\n', Nsuccess/Nloop, MED);
end