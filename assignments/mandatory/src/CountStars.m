function numStars = CountStars(I)

% Contruct our structuring element and perfom morphological opening
disk = strel('disk', 5, 4);
open = imopen(I, disk);

% Remove the background and perform some magic adjustments
F = I - open;
F = imadjust(F, [0.03, 0.98]);

% Threshold the image (with magic)
bw = im2bw(F, 0.07);

% Find the connected components in the image
cc = bwconncomp(bw, 4);
numStars = cc.NumObjects;

end
