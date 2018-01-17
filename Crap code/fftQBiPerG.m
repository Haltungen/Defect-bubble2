function out = fftQBiPerG(r, y, alp, k, NN, N)

theta = linspace(0,2*pi,N+1);
theta = theta(1:end-1);
dtheta = 2*pi/N;
x1 = r*cos(theta);
x2 = r*sin(theta);
d = 1; % Periodicity in x- and y- directions

G = ops.GBiPeriodic(k, x1-y(1), x2-y(2), d, d, alp);

pos = fft(G);
neg = conj(fft(conj(G)));
out = dtheta/(2*pi).*[fliplr(neg(2:NN+1)), pos(1:NN+1)];
end
