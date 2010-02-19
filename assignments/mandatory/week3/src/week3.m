function week3( ~ )
%Week 3: SIP Mandatory assignment
%   Mandatory assignment for Signal Image Processing
%   May include some supplementary assignments as well.

    function G = AddBorder(I, n)
        % Add a black border to an image
        G = zeros(size(I,1) + 2*n, size(I,2) + 2*n);
        G = mat2gray(G);
        G = im2uint8(G);
        
        G(n + 1:size(I,1) + n, n + 1:size(I,2) + n) = I(:,:);
    end

    function G = ExpandImage(I, n)
        % Continues edge pixel to the specified size
        G = AddBorder(I, n);
        
        row_copy = I(1,:);
        for i = 1:n
            for j = 1:length(row_copy)
                G(i, n + j) = row_copy(j);
            end
        end
        
        row_copy = I(size(I,1),:);
        for i = size(I,1) + 1 + n:size(I,1) + 2*n
            for j = 1:length(row_copy)
                G(i, n + j) = row_copy(j);
            end
        end
        
        col_copy = G(:,n + 1);
        for i = 1:n
            for j = 1:length(col_copy)
                G(j, i) = col_copy(j);
            end
        end
        
        col_copy = G(:,size(G,2) - n);
        for i = size(G,2) - n:size(G,2)
            for j = 1:length(col_copy)
                G(j, i) = col_copy(j);
            end
        end
    end

    function G = MyConv2(I, filter)
        % Assume that the image have been padded in some way
        % as we cut a piece off
        
        % Get image and filter sizes
        [iM iN] = size(I);
        [fM fN] = size(filter);
        
        % init size for output image
        gM = iM - fM + 1;
        gN = iN - fN + 1;
        G = zeros(gM, gN);
        
        % Actual convolution
        % For every pixel in output image
        for gm = 1:gM
            for gn = 1:gN
                % Init val to zero
                val = 0;
                % Walk through the support and multiply the original
                % pixel by the corresponding filter value
                for fm = 1:fM
                    for fn = 1:fN
                        val = val + filter(fm, fn) * I(gm - fm + fM, gn - fn + fN);
                    end
                end
                G(gm, gn) = val;
            end
        end
    end

    function G = BoxFilter(I, n)
        N = 2*n + 1;
        M = N;
        filter = (1/(M*N)) * ones(M, N);
        
        G = ExpandImage(I, n);
        
        %G = conv2(G, filter, 'valid');
        G = MyConv2(G, filter);
        G = mat2gray(G);
        G = im2uint8(G);
        
    end

    function MultiplyPlots( ~ )
        % Stupid plots
        n = 1:10;
        figure, hold on;
        plot(n, 2*(512^2).*(2.^n));
        plot(n, 6*(512^2)*log2(512)+512^2, 'r-');
        hold off;
    end

    function G = FilterImage(I, d0, band, method)
       
        % Create filter
        if strcmp(method, 'ideal')
            h = MakeIdealFilter(I, d0);
        elseif strcmp(method, 'butter')
            n = 2;
            h = MakeButterworthFilter(I, d0, n);
        elseif strcmp(method, 'emphasis')
            alpha = 2;
            beta = 70;
            h = MakeHighFrequencyEmphasisFilter(I, d0, alpha, beta);
        else
            % Fail is fail
            error('No such method');
        end
        
        if strcmp(band, 'low')
            % Do nothing
        elseif strcmp(band, 'high')
            h = (1 - h);
        else
            % Fail
            error('Select valid band')
        end
        
        % Inspect the filter
        %figure, imshow(h, []);
        
        % Save the filter to disk
        imwrite(h, ['../report/images/' method '_' band '_filter_' num2str(d0) '.png'], 'png');
        
        % Transform image
        F = fft2(I);
        F = fftshift(F);
        %ftshow(F);
        
        % Pointwise multiplication with the filter        
        G = h.*F;
        
        % Inspect result (Fourier)
        %ftshow(G);
        
        % Reverse transform to spatial
        G = abs(ifft2(fftshift(G)));
    end

    function h = MakeHighFrequencyEmphasisFilter(I, d0, alpha, beta)
        [im in] = size(I);
        
        h = zeros(size(I));
        
        for u = 1:im
            for v = 1:in
                D = sqrt((u - im/2)^2 + (v - in/2)^2);
                if D > 0
                    h(u, v) = alpha + beta*(1/(1 + (D/d0)^(2)));
                elseif D == 0;
                    h(u, v) = 1;
                end
            end
        end   
    end

    function h = MakeIdealFilter(I, d0)
        [im in] = size(I);
        
        h = zeros(size(I));
        
        for u = 1:im
            for v = 1:in
                D = sqrt((u - im/2)^2 + (v - in/2)^2);
                if D < d0
                    h(u, v) = 1;
                end
            end
        end
    end

    function h = MakeButterworthFilter(I, d0, n)        
        [im in] = size(I);
        
        h = zeros(size(I));
        
        for u = 1:im
            for v = 1:in
                D = sqrt((u - im/2)^2 + (v - in/2)^2);
                h(u, v) = 1/(1 + (D/d0)^(2*n));
            end
        end
    end

    function ftshow(F)
       figure, imshow(log(abs(F)), []); 
    end

    function Repair(I)
       F = fft2(I);
       F = fftshift(F);
       
       ftshow(F);
       
       % Remove noise in noisy.tiff
       %F(55:70, 160:180) = 0;
       %F(180:200, 70:90) = 0;
       
       % Remove noise in berlinger.tiff
       F(370:395, 10:70) = 0;
       F(370:395, 260:290) = 0;
       F(120:135, 220:260) = 0;
       F(115:135, 470:500) = 0;      

       
       ftshow(F);
       G = abs(ifft2(fftshift(F)));
       
       % Repair berlinger.tiff
       G = FilterImage(G, 100, 'low', 'butter');
       
       figure, imshow(G, []);
       
    end

    function C = CostSpatial(M, N)
        C = 2*M^2 * N;
    end

    function C = CostFrequency(M)
        C = 6*(M^2)*log2(M)+(M^2);
    end

    function C = CostFrequencyA(M)
        C = 4*(M^2)*log2(M)+(M^2);
    end

    function FindCosts(M)
        for N = [4 8 16 32]
            s = CostSpatial(M, N)
            f = CostFrequency(M)
        end
    end

    function MinFilterSize( ~ )
        max = 12;
        miD = zeros(1, max);
        miC = zeros(1, max);
        Ms = zeros(1, max);
        for m = 1:max
            Ms(m) = 2^m;
            miD(m) = 2^(floor(log2(3*log2(2^m) + 0.5)));
            miC(m) = 3*m + 0.5;
        end
        close all;
        figure;
        stairs(Ms, miD);
        figure;
        plot(Ms, miC);
    end

    function run( ~ )
        %g1 = imread('../../../../images/lenna.tiff');
        %g1 = imread('../../../../images/noisy.tiff');
        %g1 = imread('../../../../images/berlinger.tiff');
        %g1 = imread('../../../../images/square.tiff');
        g1 = imread('../../../../images/unix.tiff');
        %imwrite(g1, '../report/images/unix.png', 'png');
        imshow(g1, []);
        
        %Repair(g1);
        
        %FindCosts(1024);
        
        %MinFilterSize();
        
        G1 = FilterImage(g1, 45, 'low', 'butter');
        %G2 = FilterImage(g1, 45, 'high', 'ideal');
        %imwrite(G1, [gray], '../report/images/unix_butter_result_45.png', 'png');
        %imwrite(G2, [gray], '../report/images/unix_high_result_45.png', 'png');
        figure, imshow(G1, []);
        %figure, imshow(G2, []);
    end

run();

end