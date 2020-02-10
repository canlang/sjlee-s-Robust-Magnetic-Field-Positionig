% A = imread('N1-7F.png');
% xWorldLimits = [0 1670/20];
% yWorldLimits = [0 660/20];
% RA = imref2d(size(A),xWorldLimits,yWorldLimits);
% figure
% imshow(A,RA);
% zero point (20, 659), 1/20 scaled image

A = imread('map/N1-7F.png','BackgroundColor',[1 1 1]);

xWorldLimits = [-1 1650/20];
yWorldLimits = [-1 660/20];
RA = imref2d(size(A),xWorldLimits,yWorldLimits);
imshow(flipud(A),RA);
axis xy;

hold on
plot(est(:,1), est(:,2),'*r')