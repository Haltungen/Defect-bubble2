function out = f_hat(func, r, NN, N)

theta = linspace(0,2*pi,N+1);
theta = theta(1:end-1);
dtheta = 2*pi/N;
x1 = r*cos(theta);
x2 = r*sin(theta);

func_vals = feval(func,x1,x2);

pos = fft(func_vals);
neg = conj(fft(conj(func_vals)));
out = dtheta/(2*pi).*[fliplr(neg(2:NN+1)), pos(1:NN+1)];
end
