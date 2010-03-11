function I = ToIm(I)

% Make sure that matlab see this as an image
I = mat2gray(I);
I = im2uint8(I);

end
