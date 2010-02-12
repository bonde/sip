function week2( ~ )
%Week 2: SIP Mandatory assignment
%   Mandatory assignment for Signal Image Processing
%   May include some supplementary assignments as well.

    function W = myW(n, N)
       W = exp((2*pi*1i*(n))/N);
    end

    function G = IterativeFFT(g)
        % This one fucks up --- it's borked
        % It's from Cormen et al.
        G = bitrevorder(g);
        N = length(g);
        for s = 1:log(N)
            M = 2^(s);
            Wm = exp(2*pi*1i/M);
            for k = 1:M:N
                W = 1;
                for j = 0:M/2 - 1
                    u = G(k + j);
                    wS = W*G(k + j + M/2);

                    G(k + j) = u + wS;
                    G(k + j + M/2) = u - wS;
                    W = W*Wm;
                end
            end
        end
    end

    function F = IterativeFFT2(g)
        [M N] = size(g);
        F = double(g(:,:));

        for u = 1:M
            F(u,1:N) = IterativeFFT(F(u,1:N));
        end
        %for v = 1:N
        %    F(1:M,v) = IterativeFFT(F(1:M,v));
        %end
    end

    function G = RecursiveFFT(g, N)
        % Recursive divide and conquer algorithm from CLRS p. 835
        %
        % O(Nlog(N)) but with very high constants because of
        % list allocation.
        if N == 1
            % Trivial
            G = g;
            return
        else
            G  = zeros(size(1:N));
            ge = zeros(size(1:(N/2)));
            go = zeros(size(1:(N/2)));
            
            % Divide
            for v = 1:(N/2)
                ge(v) = g(2*v);
                go(v) = g(2*v - 1);
            end
    
            % Conquer
            Ge = RecursiveFFT(ge, (N/2));
            Go = RecursiveFFT(go, (N/2));
            
            % Combine
            for v = 1:(N/2)
                wS = myW(-v + 1, N)*Go(v);
                G(v)         = Ge(v) + wS;
                G(v + (N/2)) = Ge(v) - wS;
            end
        end
    end

    function F = RecursiveFFT2(I)
        % 2-D recursive Fourier transform
        [M N] = size(I);
        F = double(I(:,:));

        for u = 1:M
            F(u,1:N) = RecursiveFFT(F(u,1:N), N);
        end
        for v = 1:N
            F(1:M,v) = RecursiveFFT(F(1:M,v), M);
        end
    end

    function run( ~ )
        g1 = imread('../../../../images/square.tiff');
        %g = imread('../../../../images/lenna.tiff');
        g2 = imread('../../../../images/noisy2.tiff');
        
        %FgR = RecursiveFFT2(g);
        %FgR = fftshift(FgR);
        
        %FgI = IterativeFFT2(g);
        %FgI = fftshift(FgI);
        
        figure, imshow(g1);
        figure, imshow(g2);
        %figure, imshow(log(abs(1 + FgI)), []);
        %figure, imshow(log(abs(1 + FgR)), []);
        
        % Control
        ff1 = fft2(g1);
        ff1 = fftshift(ff1);
        ff1 = log2(abs(1 + ff1));
        figure, imshow(ff1, []);
        ff2 = fft2(g2);
        ff2 = fftshift(ff2);
        ff2 = log2(abs(1 + ff2));
        figure, imshow(ff2, []);
        
        %imwrite(g1, '../report/images/img1.png', 'png');
        %imwrite(g2, '../report/images/img2.png', 'png');
        
        % Does not work properly, damn you matlab
        %imwrite(ff1, [hot], '../report/images/fft1.png', 'png');
        %imwrite(ff2, [hot], '../report/images/fft2.png', 'png');
    end

run();

end

