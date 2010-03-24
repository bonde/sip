function I = Reduce(I, n)
% Assume 255 color, i.e. the image colors have been adjusted
% Basic idea:
%       Like a histogram, put the pixel values in
%       bins according to their value.
%       Each bin will have its own color.

% Making bins the lol-way
int_width = 255/n;
bins(1,1) = 0;
bins(1,2) = int_width;
if n > 1
    for i = (2:n)
        bins(i, 1) = bins(i - 1, 2);
        bins(i, 2) = bins(i, 1) + int_width;
    end
end

% Rambo-fix for ensuring that max values are put in the
% last bin
bins(n, 2) = bins(n, 2) + 1;

% Get the new color values for the bins
vals = floor(mean(bins, 2));

% Put pixel values in bins, thus assigning new colors
for x = (1:size(I,1))
    for y = (1:size(I,2))
        for i = 1:n
            if bins(i,1) <= I(x,y) && I(x,y) < bins(i,2)
                I(x,y) = vals(i);
                break;
            end
        end
    end
end

end
