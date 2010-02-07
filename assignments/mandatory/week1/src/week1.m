function week1( ~ )
%Week 1: SIP Mandatory assignment
%   Mandatory assignment for Signal Image Processing
%   Includes some supplementary assignments as well.

    function M = adjust(I)
        % Supplementary assignment 1.3

        % Rather use a confidence interval
        % The mean would have been acceptable also
        confidence = [0.05; 0.95];
        img_min = quantile(quantile(I, confidence(1)), confidence(1));
        img_max = quantile(quantile(I, confidence(2)), confidence(2));

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

    function I = reduceVar(I, n)
        % Uses bins from the imhist()-function
        % but does the same thing as reduce( ... )
        [count, bins] = imhist(I, n + 1);
        bins(n + 1) = bins(n + 1) + 1;
        
        % Put pixel values in bins, thus assigning new colors
        for m = 1:size(I,1)
            for n = 1:size(I,2)
                for i = 1:(length(bins) - 1)
                    if bins(i) <= I(m,n) && I(m,n) < bins(i + 1)
                        I(m,n) = floor((bins(i) + bins(i + 1))/2);
                        break;
                    end
                end
            end
        end
    end

    function I = factorResize(I, factor)
        % Resize image with a 1/factor scale
        % Each pixel is given the mean of the factor*factor
        % neighbourhood
        % -------------
        % | x |   |   |
        % |-----------|
        % |   |   |   |    <- Mean on this given to x
        % |-----------|
        % |   |   |   |
        % |-----------|
        %
        
        [m, n] = size(I);
        
        % Calculate the extra width required for no modulo
        extra_width = factor - mod(n, factor);
        extra_height = factor - mod(m, factor);
        
        if (extra_width/factor) < 0.5 && mod(n, factor) ~= 0
            tmp = mat2gray(zeros(m, n + extra_width));
            tmp = im2uint8(tmp);
            tmp(1:m, 1:n) = I(1:m,1:n);
            
            cop_col = I(1:m, n);
            for i = 1:length(cop_col)
                for j = n + 1:n + extra_width
                    tmp(i, j) = cop_col(i);
                end
            end
            
            % Set the original image the modified
            I = tmp;
            
            % Get rid of the extra image
            clearvars tmp;
        elseif extra_width ~= 0
            % We want to clip the image
            I = I(1:m , 1:n - mod(n, factor));
        end
        
        [m, n] = size(I);
        
        if (extra_height/factor) < 0.5 && mod(m, factor) ~= 0
            tmp = mat2gray(zeros(m + extra_height, n));
            tmp = im2uint8(tmp);
            tmp(1:m, 1:n) = I(:,:);
            
            cop_row = I(m, 1:n);
            for i = m + 1:m + extra_height
                for j = 1:length(cop_row)
                    tmp(i, j) = cop_row(j);
                end
            end
            I = tmp;
            clearvars tmp;
        elseif extra_height ~= 0
            I = I(1:m - mod(m, factor), 1:n);
        end
        
        [m, n] = size(I);

        for x = 1:factor:n
            for y = 1:factor:m
                I(y:y + factor - 1, x:x + factor - 1) = floor(mean2(I(y:y + factor - 1, x:x + factor - 1)));
            end
        end

    end

    function I = quatersize(I, times)
        % Supplementary excercise 1.4
        for x = 1:times
            I = I(1:2:size(I,1), 1:2:size(I,2));
        end
    end

    function I = quadouble(I, times)
        %Supplementary excercise 1.4
        for i = (1:times)
            M = zeros(size(I) * 2);
            xoff = 0;
            
            for x = (1:size(I, 1))
                yoff = 0;
                for y = (1:size(I, 2))
                    M(x + xoff, y + yoff) = I(x,y);
                    M(x + 1 + xoff, y + yoff) = I(x,y);
                    M(x + 1 + xoff, y + 1 + yoff) = I(x,y);
                    M(x + xoff, y + 1 + yoff) = I(x,y);
                    yoff = yoff + 1;
                end
                xoff = xoff + 1;
            end
            M = mat2gray(M);
            M = im2uint8(M);
            I = M;  
        end
    end

    function I = resolution(I, scale)
        kernel = fspecial('gaussian', [3 3], 5);
        I = imfilter(I, kernel, 'symmetric', 'conv');
        I = I(1:2:size(I,1), 1:2:size(I,2));
    end
        

    function run( ~ )
        filepath = '../../../../images/cmp1.gif';
        %filepath = '../../../../images/idotyl.tiff';

        I = imread(filepath);

        I = adjust(I);

        n = 8;
        %R = reduce(I, n);
        R = reduceVar(I, n);
        
        factor = 2;
        M = factorResize(R, factor);
        %M = I(1:factor:size(M,1), 1:factor:size(M,2));
        
        imshow(I), figure, imhist(I, 128);
        figure, imshow(M), figure, imhist(M, n);
        %figure, imshow(R);
    end

run();

end

