
Nloop = 2;
Nfailure = 0;
Nsuccess = 0;
errs = cell(Nloop,1);
for i=1:Nloop
    close all
    err = ILoA('KI-1F','MATE20pro',3,.8,1);
    %%    
    if sum(err<3)>length(err)/2
        disp converged!
        Nsuccess = Nsuccess+1;
        errs{i} = err;
    else
        Nfailure = Nfailure+1;
        disp fail!
    end
end