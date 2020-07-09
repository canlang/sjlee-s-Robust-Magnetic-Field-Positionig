close all;
n_particle = 4000;
site = 'KI-1F';          % candidates: N1-7F (iphone), KI-1F (iphone, MATE20pro, S9)
device = 'S9';
trj = 2;
visualize = true;
l = 1;  
intp_level = .8;
err = ILoA(n_particle,site,device,trj,intp_level,visualize,l);                