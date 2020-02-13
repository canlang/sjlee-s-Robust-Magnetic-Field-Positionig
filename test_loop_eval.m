
Nloop = 1;

% persistent intp;
% intp = [.1, .2, .3, .5, .8, 1.0, 1.2];
% intp = [.2, .3, .5, .8, 1.0, 1.2];
intp = [.8];
errs = cell(Nloop,length(intp));
trj_idx = 3;
device_name = 'S9';
site_name = 'KI-1F';

for j=1:length(intp)
    Nfailure = 0;
    Nsuccess = 0;
    Dintp = intp(j);
    for i=1:Nloop
        close all
        err = ILoA('KI-1F','S9',trj_idx,Dintp,true);
        %%    
        if sum(err<5)>length(err)/2
            % disp converged!
            Nsuccess = Nsuccess+1;
            errs{j,i} = err;
        else
            Nfailure = Nfailure+1;
            % disp fail!
        end
    end
    all_errs = errs{j,:};
    MED = mean(rmoutliers(all_errs));
    fprintf('GeoMapInt: %.1f, Convergence rate: %.1f, MED: %.2f\n'...
    , intp(j), Nsuccess/Nloop*100, MED);
    % return
end
% save(sprintf('est-result/%s-s%d-%s-errs.mat',site_name,trj_idx,device_name),'errs')
%%
% isempty(