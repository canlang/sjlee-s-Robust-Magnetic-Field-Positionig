% x = -pi:0.01:pi;
% % conv(sin(x),cos(x))
% % plot(conv(sin(x),cos(x)))
% 
% 
% [C1,lag1] = xcorr(cos(x),sin(x)); 
% % figure
% % ax(1) = subplot(2,1,1); 
% plot(0,C1,'k')
% ylabel('Amplitude')

clc,clf, clear, close all

load('mats/rotated_magnetic_field.mat')

std_C = stdfilt(C(1,:));
[~,maxI] = max(std_C);


title_idx = 'xyz';
figure
for i=1:3
    x = C(i,1:90);
    y = C(i,100:end);

    subplot(4,1,i)
    plot(x)
    hold on
    plot(fliplr(y))
    title(sprintf('trends of %s magnetic comp. ',title_idx(i)));
    hold off
end
 
x = C(1,1:90);
y = C(1,100:end);

% R_xy = ifft(fft(x,length(x)) * conj(fft(y,length(y))))
% [C1,lag1] = xcorr(x,y);
subplot(414)
plot(zscore(x))
hold on 
plot(fliplr(zscore(y)))
hold off
% plot(lag1, C1);

set(gcf,'units','points','position',[500,500,1000,600])
sdf(gcf,'sj2')
%%
fun = @(X) dtw(x,fliplr(y)+X);
opt_arg = fminsearch(fun,0);
figure
dtw(x,fliplr(y)+opt_arg)
set(gcf,'units','points','position',[500,500,1000,600])
sdf(gcf,'sj2')