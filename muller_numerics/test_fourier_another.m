% Illustrates how the singularity of the normal derivative of the
% quasi-periodic greens function depends on the direction pf approach

clear all
close all

theta = 0:0.1:2*pi;
M = length(theta);
d = 1;
NN = 5;
m = 5;
n = (-NN:NN).';
N = 5;
Ntilde = 100;
R_b = 0.5;
r_p = 0.48;
theta_p = theta(floor(end/2))+0.1;
y = r_p*[cos(theta_p), sin(theta_p)];
alpha = [pi/2, pi/8];
k1 = 1;
k2 = 0.8;
figure
hold on

for r = 0.48:0.01:0.52
y = r*[cos(theta_p),sin(theta_p)];
x1 = R_b*cos(theta);
x2 = R_b*sin(theta);
G_ewald = tools.GBiPeriodic(k1, x1-y(1), x2-y(2), d, d, alpha);%-ops.GBiPeriodic(k2, x1-y(1), x2-y(2), d, d, alpha);
G_ewald_nu = dot([cos(theta);sin(theta)],tools.GradGBiPeriodic(k1, x1-y(1), x2-y(2), d, d, alpha));
G_0 = -1i/4*besselh(0,1,k1*sqrt((x1-y(1)).^2+(x2-y(2)).^2));
rprime = sqrt((x1-y(1)).^2+(x2-y(2)).^2);
G_0_nu = -1i/4*k1*0.5*(besselh(-1,1,k1*rprime)-besselh(1,1,k1*rprime)).*dot([x1-y(1);x2-y(2)],[cos(theta);sin(theta)])./rprime;
plot(theta,real(G_ewald_nu))%,theta,imag(G_ewald_nu));
end
hold off