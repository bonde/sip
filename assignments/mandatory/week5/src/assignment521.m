function assignment521 ( ~ )

close all;

% Images
%[g1 cmap] = imread('../../../../images/lenna.tiff', 'tiff');
[g1 cmap] = imread('../../../../images/R1.tiff');
%g1 = imread('../../../../images/noisy.tiff');
%g1 = imread('../../../../images/berlinger.tiff');
%[g1 cmap] = imread('../../../../images/square.tiff');
%g1 = imread('../../../../images/unix.tiff');

S = SignalWeek5();

padding = 1000;
GS = ExpandSignal(S, padding);
GS = Gauss(GS, 5/4);
N = length(S);

Dx = ImDerivative(GS, 'dx','b', 1);
Dx = Dx(:, padding+1:N+padding);
Dx2 = ImDerivative(GS, 'dx','f', 2);
Dx2 = Dx2(:, padding+1:N+padding);
[C locs] = ZeroCrossings(Dx2);

hold on;
figure(1);
plot([1:N] - N, zeros([1, N]), 'color', [0.6 0.6 0.6]);
plot([1:N]-N, Dx, 'g');
plot([1:N] - N, Dx2, 'm');
plot(locs - N, zeros(size(locs)), 'xr');
ylim([-1 1]);
hold off;

GS = GS(:, padding+1:N+padding);
[vals locs] = EdgeDetect(GS, 0.8);

figure(2);
plot([1:N]-N, S, 'color', [1 0.4 0.7]);
%plot([1:N]-N, S);
hold on;
plot([1:N]-N, GS);
plot(locs-N, vals, 'or');
xlabel('Depth in meters')
ylabel('Temperature in Degree Celsius')
hold off;

end
