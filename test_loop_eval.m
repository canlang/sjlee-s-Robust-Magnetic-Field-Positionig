
Nloop = 100;
% persistent intp;
% intp = [.2, .3, .5, .8, 1.0, 1.2];
% intp = [.2, .3, .5, .8, 1.0, 1.2];
intp = [.8];

for k=1:3       % There is three trajectory (scenarios) for testing
    errs = cell(length(intp),Nloop);
    trj_idx = k;
    device_name = 'S9';
    site_name = 'KI-1F';        % candidates: N1-7F (iphone), KI-1F (iphone, MATE20pro, S9)

    for j=1:length(intp)
        Nfailure = 0;
        Nsuccess = 0;
        Dintp = intp(j);
        parfor i=1:Nloop           % available 'parfor' for speed up
            close all
            err = ILoA(site_name,device_name,trj_idx,Dintp,false);
            %%    
            if sum(err<5)>length(err)/2
                % disp converged!
                Nsuccess = Nsuccess+1;
                errs{j,i} = err;
            else
                Nfailure = Nfailure+1;
                errs{j,i} = err;
                % disp fail!
            end
        end
        all_errs = errs{j,:};
        MED = mean(rmoutliers(all_errs));
        fprintf('GeoMapInt: %.1f, Convergence rate: %.1f, MED: %.2f\n'...
        , intp(j), Nsuccess/Nloop*100, MED);
        % return
    end
    save(sprintf('est-result/%s-s%d-%s-errs-rep(%d)-intp(%.1f).mat',site_name,trj_idx,device_name,Nloop,intp),'errs')
end
%%
