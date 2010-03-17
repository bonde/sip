function G = AddBorder(I, n)
% Add a black border to an image
G = zeros(size(I,1) + 2*n, size(I,2) + 2*n);
G = mat2gray(G);
G = im2uint8(G);

G(n + 1:size(I,1) + n, n + 1:size(I,2) + n) = I(:,:);

end
