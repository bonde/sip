function week4( ~ )
%Week 3: SIP Mandatory assignment
%   Mandatory assignment for Signal Image Processing
%   May include some supplementary assignments as well.

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

    function ftshow(F)
       figure, imshow(log(abs(F)), []); 
    end

% PART 1

    function LP = AlternativeLaplacianPyramid(I, P, sigma)
       % Construct the Laplacian pyramid differently
       GP = GaussianPyramid(I, P, sigma);
       LP = cell(P, 1);

       for p = 1:P - 1
           LP{p, 1} = GP{p, 1} - Gauss(GP{p, 1}, sigma);
       end

       LP{P, 1} = GP{P, 1};
    end

    function LP = LaplacianPyramid(I, P, sigma)
        % Generate the Laplacian pyramid to a certain level
        GP = GaussianPyramid(I, P, sigma);

        LP = cell(P, 1);

        for p = 1:P - 1
            LP{p, 1} = GP{p, 1} - Upsample(GP{p + 1, 1});
        end

        LP{P, 1} = GP{P, 1};
    end

    function G = Upsample(I)
       % Upsamples an image with factor two in the Fourier domain
       F = fftshift(fft2(I));
       offM = size(I, 1)/2;
       offN = size(I, 2)/2;

       G = zeros(size(I)*2);
       G(offM + 1:size(G,1) - offM, offN + 1:size(G,2) - offN) = F(:,:);
       G = abs(ifft2(G));
       G = double(G);
    end

    function G = Downsample(I)
       % We expect that the image have been blurred
       F = fftshift(fft2(I));

       % Check if we are one-dimensional
       if size(I, 1) ~= 1
           offM = size(I,1)/4;
       else
           offM = 0;
       end

       if size(I, 2) ~= 1
           offN = size(I,2)/4;
       else
           offN = 0;
       end

       % Lols, this work too, dumbass
       F = F(offM + 1:3*offM, offN + 1:3*offN);
       %F = F(offM + 1:size(I,1) - offM, offN + 1:size(I,2) - offN);

       G = abs(ifft2(F));
    end

    function GP = ReconstructFromLaplacian(LP)
        % Reconstruct the original image from a Laplacian pyramid
        % We do this by reconstructing the Gaussian Pyramid
        P = size(LP,1);
        GP = cell(P, 1);

        % The top of the pyramids is equal
        GP{P, 1} = LP{P, 1};

        for p = P:-1:2
            GP{p - 1, 1} = LP{p - 1, 1} + Upsample(GP{p, 1});
        end
    end

    function GP = GaussianPyramid(I, Q, sigma)
        % Construct a gaussian pyramid of height Q
        GP = cell(Q, 1);

        GP{1, 1} = double(I);

        for q = 1:Q - 1
            GP{q + 1, 1} = Downsample(Gauss(GP{q, 1}, sigma));
        end
    end

    function G = Gauss(I, sigma)
        % Compute the gaussian with specified sigma
        % Convolution is done in the frequency domain
        F = fft2(I);

        h = fspecial('gaussian', size(I), sigma);
        h = fft2(h);

        G = abs(fftshift(ifft2(h.*F)));
    end

    function ShowPyramid(PYR, cmap)
        figure;
        hold on;
        for i = 1:length(PYR)
            subplot(2, 3, i), imshow(PYR{i, 1}, cmap);
        end
        hold off;
    end

    function part1( ~ )
        [g1 cmap] = imread('../../../../images/lenna.tiff', 'tiff');
        %[g1 cmap] = imread('../../../../images/R1.tiff');
        %g1 = imread('../../../../images/noisy.tiff');
        %g1 = imread('../../../../images/berlinger.tiff');
        %g1 = imread('../../../../images/square.tiff');
        %g1 = imread('../../../../images/unix.tiff');

        sigma = 4;
        levels = 6;
        % Write something
        i = Downsample(g1);
        figure, imshow(i, [])
    end

% PART 2

    % Constructs the random signal
    function S = ConstructSignal( ~ )
        N = 64;
        S = [zeros(1,20),ones(1,3),zeros(1,2),ones(1,5)];
        S = [S,zeros(1,N-length(S))] + 0.1*randn(1,N);
    end

    % Scale the signal to a certain amount of levels
    function I = ScalesOfSignal(S, levels)
        I = zeros(levels, length(S));
        I(1,:) = S;
        for sigma = 1:levels - 1
            I(sigma + 1, :) = Gauss(S, sigma);
        end
    end

    % Same as above, but in steps of exponentials of two
    function I = AltScalesOfSignal(S, levels)
        I = zeros(levels, length(S));
        I(1,:) = S;
        for sigma = 1:levels - 1
            I(sigma + 1, :) = Gauss(S, 2^(sigma - 1));
        end
    end

    % "Pixelwise" differerence
    function e = DetectDiff(S)
        e = zeros(size(S));
        for n = 1:size(S,2) - 1
            e(n) = S(n) - S(n + 1);
        end

    end

    function C = CountPositive(S)
        C = 0;
        for i = 1:length(S)
            if sign(S(i)) > 0
                C = C + 1;
            end
        end
    end

    % Remove zero values in a pair of vals and locations
    function [nVals nLocs] = PruneZeroes(vals, locs)
        nVals = zeros(CountPositive(vals), 1);
        nLocs = zeros(size(nVals));

        j = 1;
        for i = 1:length(vals)
            if sign(vals(i)) == 1
                nVals(j) = vals(i);
                nLocs(j) = locs(i);
                j = j + 1;
            end
        end

    end

    % Returns the number and location of zero crossings in a signal
    function [C locs] = ZeroCrossings(S)
        prevsign = 0;
        C = 0;

        for n = 1:length(S)
            thissign = sign(S(n));
            if (thissign * prevsign) == -1
                C = C + 1;
                locs(C) = n;
            end
            if thissign ~= 0
                prevsign = thissign;
            end
        end
    end

    % Main
    function part2( ~ )
        S = ConstructSignal();
        S = adjust(S);

        s1 = ScalesOfSignal(S, length(S));
        %s2 = AltScalesOfSignal(S, length(S));

        figure('Position', [1 1 600 800]);
        for i = 1:6
            FDev = DetectDiff(s1(i, :));
            SDev = DetectDiff(FDev);

            [C locs] = ZeroCrossings(FDev);
            vals = zeros(size(locs));

            subplot(3, 2, i);
            hold on;
            for j = 1:length(locs)
                tmp = s1(i, :);
                next = s1(i + 1, :);
                if sign(FDev(locs(j))) == 1
                    vals(j) = tmp(locs(j));
                end
            end

            [vals locs] = PruneZeroes(vals, locs);

            subplot(3, 2, i);
            plot(s1(i, :));
            plot(locs,vals,'x', 'Color', 'red');
            axis([0 70 0 260])

            hold off;
        end

        figure('Position', [1 1 600 800]);
        for i = 1:6

            FDev = DetectDiff(s1(i, :));
            SDev = DetectDiff(FDev);

            [C locs] = ZeroCrossings(SDev);
            vals = zeros(size(locs));

            subplot(3, 2, i);
            hold on;
            % Plot second derivate
            %plot(SDev,'Color', [1 0.7 0.7]);
            plot(FDev);
            for j = 1:length(locs)
                vals(j) = FDev(locs(j));
            end
            %scatter(locs,vals,'x');
            xlim([0 70])

            hold off;
        end

        %figure, imshow(s1, []);
        %figure, imshow(s2, []);
    end

close all;
part2();

end
