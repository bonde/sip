function M = Adjust(I, varargin)
% Supplementary assignment 1.3

% Rather use a confidence interval
% The mean would have been acceptable also
confidence = [0.05; 0.95];
img_min = quantile(quantile(I, confidence(1)), confidence(1));
img_max = quantile(quantile(I, confidence(2)), confidence(2));

% Target interval
target_min = 0;
target_max = 255;

% Init the scale
scale = (target_max - target_min) / double((img_max - img_min));

% Init target image
M = zeros(size(I));

% Very clever notation instead of a for-loop
% For every pixel in original we subtract the image minimum, then
% multiply by the calculated scale. Finally add our interval min.
M(:,:) = ((I(:,:) - img_min) * scale + target_min);

% Make sure that matlab see this as an image
M = ToIm(M);

end
