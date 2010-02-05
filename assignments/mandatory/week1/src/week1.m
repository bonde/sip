function week1( ~ )
%Week 1: SIP Mandatory assignment
%   Mandatory assignment for Signal Image Processing
%   Includes some supplementary assignments as well.

    function M = adjust(I)
        % Supplementary assignment 1.3

        % Rather use a confidence interval
        % The mean would be acceptable also
        confidence = [0.05; 0.95];
        img_min = quantile(quantile(I, confidence(1)), confidence(1));
        img_max = quantile(quantile(I, confidence(2)), confidence(2));
        %img_min = 0;
        %img_max = 120;

        % Target interval
        my_min = 0;
        my_max = 255;

        % Init the scale
        scale = (my_max - my_min) / double((img_max - img_min));

        % Init target image
        M = zeros(size(I));

        % Very clever notation instead of a for-loop
        % For every pixel in original we subtract the image minimum, then
        % multiply by the calculated scale. Finally add our interval min.
        M(:,:) = ((I(:,:) - img_min) * scale + my_min);

        % Make sure that matlab see this as an image
        M = mat2gray(M);
        M = im2uint8(M);
    end

    function I = reduce(I, n)
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
        
        % Rambo-fix
        bins(n, 2) = bins(n, 2) + 1;
        
        % Get the new color values for the bins
        vals = floor(mean(bins, 2));
        
        % Put pixel values in bins, thus assigning new colors
        for x = (1:size(I,1))
            for y = (1:size(I,2))
                for i = (1:n)
                    if bins(i,1) <= I(x,y) && I(x,y) < bins(i,2)
                        I(x,y) = vals(i);
                        break;
                    end
                end
            end
        end
    end

%filepath = '../../../images/cmp1.gif';
filepath = '../../../images/idotyl.tiff';

I = imread(filepath);

I = adjust(I);

n = 4;
I = reduce(I, n);
imshow(I);
figure, imhist(I, n)

end

