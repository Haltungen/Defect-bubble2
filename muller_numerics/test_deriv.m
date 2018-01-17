clear all
%close all

rho = 5000;
rho_b = 1;
kappa = 5000;
kappa_b = 1;
omega = 0.2;

delta = rho_b/rho;
v = sqrt(rho/kappa);
v_b = sqrt(rho_b/kappa_b);

k0 = omega * v;
kb = omega * v_b;

nPoints = 1000;
xMax = 0.07;
nX = 2^6;
hX = 2*xMax/(nX-1);
x1 = -xMax:hX:xMax;
x2 = zeros(1,nX);
alpha = [pi/2, pi/8];

R_b = 0.05;  % Size of crystal bubbles
phi = rand(nPoints,1);  
Sphi = zeros(nX,1);

for m=1:nX
    Sphi(m) = ops.S_H(k0, R_b, phi, x1(m), x2(m));
end
plot(x1,Sphi)

h = 10^-4;
h0 = 0%10^-8;
numder_plus = (ops.S_H(k0, R_b, phi, R_b+h0+h, 0)-ops.S_H(k0, R_b, phi, R_b+h0, 0))/h;
numder_minus = (-ops.S_H(k0, R_b, phi, R_b-h-h0, 0)+ops.S_H(k0, R_b, phi, R_b-h0, 0))/h;
Kphi = ops.Kstar_H_eval(k0, R_b, phi, R_b, 0);
der_plus = phi(1)/2 + Kphi;
der_minus = -phi(1)/2 + Kphi;
error_plus = abs((numder_plus-der_plus)/der_plus);
error_minus = abs((numder_minus-der_minus)/der_minus);

fprintf('Error from outside: %.5e\n', error_plus);
fprintf('Error from inside: %.5e\n', error_minus);
