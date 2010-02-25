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
       
       offM = size(I,1)/4;
       offN = size(I,2)/4;
       
       F = F(offM + 1:size(I,1) - offM, offN + 1:size(I,2) - offN);
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
        
        imshow(Gauss(g1, 8) - Gauss(g1, 16), cmap);
        %figure, imhist(g1);
        
        sigma = 4;
        levels = 6;
        
        %LP1 = LaplacianPyramid(g1, levels, sigma);
        %LP2 = AlternativeLaplacianPyramid(g1, levels, sigma);
        
        %ftshow(fftshift(fft2(g1)));
        
        %GP1 = ReconstructFromLaplacian(LP1);
        %GP2 = ReconstructFromLaplacian(LP2);
        
        %ShowPyramid(GP1, cmap);
        %ShowPyramid(LP2, cmap);
        
        %reconst1 = im2uint8(mat2gray(GP1{1,1}));
        %reconst2 = im2uint8(mat2gray(GP2{1,1}));
        
        %reconst1 = adjust(reconst1);
        %reconst2 = adjust(reconst2);
        %g1 = adjust(g1);
        
        %diff1 = im2uint8(mat2gray(g1 - reconst1));
        %diff2 = im2uint8(mat2gray(g1 - reconst2));
        %for m = 1:size(g1,1)
        %    for n = 1:size(g1,2)
        %        if g1(m,n) ~= reconst1(m,n)
                    %g1(m,n)
                    %reconst1(m,n)
        %            g1(m,n) = 255;
        %        end
        %    end
        %end
        
        %figure, imshow(g1, []);
        
        %figure, imhist(reconst1);
        %figure, imhist(reconst2);
        %imwrite(g1, '../report/images/org.png', 'png');
        %imwrite(reconst1, '../report/images/reconst1.png', 'png');
        %imwrite(reconst2, '../report/images/reconst2.png', 'png');
        
        %ftshow(fftshift(fft2(reconst1)));
        %ftshow(fftshift(fft2(reconst2)));
        
        %figure, imshow(reconst1, cmap);
        %figure, imshow(reconst2, cmap);
        
        %figure, imshow(diff1, cmap);
        %figure, imshow(diff2, cmap);
    end

% PART 2

    function S = ConstructSignal( ~ )
        N = 64;
        S = [zeros(1,20),ones(1,3),zeros(1,2),ones(1,5)];
        S = [S,zeros(1,N-length(S))] + 0.1*randn(1,N);
        plot(S)
    end

    function part2( ~ )
        S = ConstructSignal();
        
        P = zeros(length(S), length(S));
        for i = 1:size(P,1)
            P(i,:) = S;
        end
        
        %figure, imshow(S, []);
        
        F = fftshift(fft2(S));
        
        figure, plot(Gauss(S, 0.2));

        %s = 0.20;
        %n = length(S);
        %x = -1/2:1/(n-1):1/2;
        %h = exp( -(x.^2)/(2*s^2));
        %h = h / sum(sum(h));
        
        figure, plot(log(abs(F)));
        %figure, plot(log(abs(h.*F)));

        %figure, plot(h);
        
        %S = abs(ifft(h.*F));
        %figure, plot(S);

    end

    function run( ~ )
        part1;
    end

close all;
run();

end