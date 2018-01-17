% Tests the fourier expansions of the quasi-periodic greens function and
% normal derivative

clear all
close all

theta = 0:0.1:2*pi;
M = length(theta);
d = 1;
NN = 2;
N1 = 1;
N2 = 1;
n = (-NN:NN).';
Ntilde = 500;
R_b = 0.05;
r_p = 0.05;
theta_p = theta(floor(end/2))+0.05;
y = r_p*[cos(theta_p), sin(theta_p)];
alpha = [pi/2, pi/8];
k1 = 1;
k2 = 0.8;

func1 = @(x1,x2) tools.GBiPeriodic(k1, x1-y(1), x2-y(2), d, d, alpha);%-tools.GBiPeriodic(k2, x1-y(1), x2-y(2), d, d, alpha);
func2 = @(x1,x2) -1i/4*besselh(0,1,k1*sqrt((x1-y(1)).^2+(x2-y(2)).^2));
func3 = @(x1,x2) dot([cos(atan2(x2,x1));sin(atan2(x2,x1))],tools.GradGBiPeriodic(k1, x1-y(1), x2-y(2), d, d, alpha));
func4 = @(x1,x2) -1i/4*k1*0.5*(besselh(-1,1,k1*sqrt((x1-y(1)).^2+(x2-y(2)).^2))-besselh(1,1,k1*sqrt((x1-y(1)).^2+(x2-y(2)).^2))).*dot([x1-y(1);x2-y(2)],[cos(atan2(x2,x1));sin(atan2(x2,x1))])./sqrt((x1-y(1)).^2+(x2-y(2)).^2);


[Ghat, Ghat_nu] = tools.fourierQBiPerG(R_b, y, alpha, k1, NN, N1, N2);%-tools.fourierQBiPerG(R_b, y, alpha, k2, NN, N, m);
Ghat_fft = tools.f_hat(func1, R_b, NN, Ntilde).';
Ghat_nu_fft = tools.f_hat(func3, R_b, NN, Ntilde).';

G_fourier = zeros(M,1);
G_nu = zeros(M,1);
G_nu_fft = zeros(M,1);
G_fourier_fft = zeros(M,1);
for l = 1:M
    G_fourier(l) = sum(Ghat.*exp(1i.*n.*theta(l)));
    G_fourier_fft(l) = sum(Ghat_fft.*exp(1i.*n.*theta(l)));
    G_nu(l) = sum(Ghat_nu.*exp(1i.*n.*theta(l)));
    G_nu_fft(l) = sum(Ghat_nu_fft.*exp(1i.*n.*theta(l)));
end
% theta1 = pi/2;
% G_test1 = sum(Ghat.*exp(1i.*n.*theta1));
% y = R_b*[cos(theta1), sin(theta1)];
% [Ghat, ~] = tools.fourierQBiPerG(R_b, y, alpha, k1, NN, N, m);
% G_test2 = sum(Ghat.*exp(1i.*n.*theta_p));

x1 = R_b*cos(theta);
x2 = R_b*sin(theta);
G_ewald = tools.GBiPeriodic(k1, x1-y(1), x2-y(2), d, d, alpha);%-ops.GBiPeriodic(k2, x1-y(1), x2-y(2), d, d, alpha);
G_ewald_nu = dot([cos(theta);sin(theta)],tools.GradGBiPeriodic(k1, x1-y(1), x2-y(2), d, d, alpha));
G_0 = -1i/4*besselh(0,1,k1*sqrt((x1-y(1)).^2+(x2-y(2)).^2));
rprime = sqrt((x1-y(1)).^2+(x2-y(2)).^2);
G_0_nu = -1i/4*k1*0.5*(besselh(-1,1,k1*rprime)-besselh(1,1,k1*rprime)).*dot([x1-y(1);x2-y(2)],[cos(theta);sin(theta)])./rprime;

% % Numerical derivative
% h = 0.01;
% R_h = R_b+h;
% x1_h = R_h*cos(theta);
% x2_h = R_h*sin(theta);
% G_ewald_0 = ops.GBiPeriodic(k1, x1-y(1), x2-y(2), d, d, alpha);
% G_ewald_h = ops.GBiPeriodic(k1, x1_h-y(1), x2_h-y(2), d, d, alpha);
% G_nu_num = (G_ewald_h-G_ewald_0)/h;

plot(n,real(Ghat_fft),n, imag(Ghat_fft),n, real(Ghat),'--*',n,imag(Ghat),'--*')
legend('real(G_{fft})', 'imag(G_{fft})','real(G_{mpl})', 'imag(G_{mpl})');
figure
plot(n,real(Ghat_nu_fft),n, imag(Ghat_nu_fft),n, real(Ghat_nu),'--*',n,imag(Ghat_nu),'--*')
legend('real(G_{\nu fft})', 'imag(G_{\nu fft})','real(G_{\nu mpl})', 'imag(G_{\nu mpl})');
figure
plot(theta,real(G_fourier),'--*',theta,imag(G_fourier),'--*',theta,real(G_ewald),theta,imag(G_ewald));
legend('real(G_{fourier})','imag(G_{fourier})','real(G_{ewald})','imag(G_{ewald})');
figure
plot(theta,real(G_nu),'--*',theta,imag(G_nu),'--*',theta,real(G_ewald_nu),theta,imag(G_ewald_nu));
% hold on
% plot(theta,real(G_nu_num),theta,imag(G_nu_num))
% hold off
legend('real(G_{\nu,fourier})','imag(G_{\nu,fourier})','real(G_{\nu,ewald})','imag(G_{\nu,ewald})','real(G_{\nu,num})','imag(G_{\nu,num})');