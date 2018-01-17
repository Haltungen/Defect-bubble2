% Not really usefull

clear all
close all

rho = 5000;
rho_b = 1;
kappa = 5000;
kappa_b = 1;
omega = 0.2;
k = omega*sqrt(rho/kappa);
k_b = omega*sqrt(rho_b/kappa_b);
delta = rho_b/rho;
nPoints = 2^6;
NN = 10;
N = 2*NN;
R_b = 0.5;
r_p = 0.4;
theta_p = 0;
alpha = [pi/2, pi/8];
y = r_p*exp(1i*theta_p);

makeF(R_b,y, k, k_b, alpha, delta, NN, N)