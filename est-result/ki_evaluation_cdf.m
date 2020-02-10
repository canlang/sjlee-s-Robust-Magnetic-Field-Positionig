close all;clearvars
load ki-err-results-at-3-sites.mat
cdfplot(a)
hold on
cdfplot(b)
cdfplot(c)
hold off
legend('s#1','s#2','s#3')