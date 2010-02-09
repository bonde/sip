function week2( ~ )
%Week 2: SIP Mandatory assignment
%   Mandatory assignment for Signal Image Processing
%   May include some supplementary assignments as well.

    function W = myW(n, N)
       W = exp((2*pi*1i*(n))/N);
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
        g = imread('../../../../images/square.tiff');
        %g = imread('../../../../images/lenna.tiff');
        %g = imread('../../../../images/noisy2.tiff');
        
        Fg = RecursiveFFT2(g);
        Fg = fftshift(Fg);
        
        imshow(g);
        figure, imshow(log(abs(1 + Fg)), []);
        
        % Control
        %ff = fft2(g);
        %ff = fftshift(ff);
        %figure, imshow(log(abs(1 + ff)), []);
    end

run();

end

