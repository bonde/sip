function numStars = CountStars(I, cmap)

disk = strel('disk', 5, 4);
open = imopen(I, disk);

F = I - open;
F = imadjust(F, [0.03, 0.98]);
bw = im2bw(F, 0.07);

figure, imshow(I, []);
figure, imshow(bw, []);

F = double(F);
figure, imshow(F, cmap);

cc = bwconncomp(bw, 4);
numStars = cc.NumObjects;

end
