function week5 ( ~ )

    function S = GenerateSignal( ~ )
        N = 64;
        S = [zeros(1,20),ones(1,3),zeros(1,2),ones(1,5)];
        S = [S,zeros(1,N-length(S))] + 0.1*randn(1,N);
        S = cumsum(S);
    end

close all;
hold off;

% Images
%[g1 cmap] = imread('../../../../images/lenna.tiff', 'tiff');
[g1 cmap] = imread('../../../../images/R1.tiff');
%g1 = imread('../../../../images/noisy.tiff');
%g1 = imread('../../../../images/berlinger.tiff');
%[g1 cmap] = imread('../../../../images/square.tiff');
%g1 = imread('../../../../images/unix.tiff');

S = GenerateSignal();
padding = 1000;
GS = ExpandSignal(S, padding);
GS = Gauss(GS, 5/4);
GS = GS(:, padding+1:N+padding);
N = length(S);

hold on;
figure(1);
plot([1:N] - N, zeros([1, N]), 'color', [0.6 0.6 0.6]);
Dx = ImDerivative(GS, 'dx','f',1);
plot([1:N]-N, Dx, 'g');
D2x = ImDerivative(Dx, 'dx','f',1);
[C locs] = ZeroCrossings(D2x);
plot([1:N] - N, D2x,'m');
plot(locs - N, zeros(size(locs)), 'xr');
ylim([-1 1]);
hold off;

figure(2);
plot([1:N]-N, S, 'color', [1 0.4 0.7]);
hold on;
plot([1:N]-N, GS);
[vals locs] = EdgeDetect(GS, 0.75);
plot(locs-N, vals, 'or');
xlabel('Depth in meters')
ylabel('Temperature in Degree Celsius')
hold off;

end
