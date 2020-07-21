%# create stacked images (I am simply repeating the same image 5 times)
img = load('clown');
I = repmat(img.X,[1 1 5]);
cmap = img.map;

A = imread('N1-7F.png','BackgroundColor',[1 1 1]);

%# coordinates
[X,Y] = meshgrid(1:size(I,2), 1:size(I,1));
Z = ones(size(I,1),size(I,2));

%# plot each slice as a texture-mapped surface (stacked along the Z-dimension)
for k=1:size(I,3)
    surface('XData',X-0.5, 'YData',Y-0.5, 'ZData',Z.*k*2, ...
        'CData',A, 'CDataMapping','direct', ...
        'EdgeColor','none', 'FaceColor','texturemap')
end
colormap(cmap)
view(3), box on, axis tight square
set(gca, 'YDir','reverse', 'ZLim',[0 size(I,3)+1])

axis image