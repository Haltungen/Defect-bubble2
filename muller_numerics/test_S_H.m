clear all
close all

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
xMax = 0.05;
nX = 20^2;
hX = 2*xMax/(nX-1);
x = -xMax:hX:xMax;

R_b = 0.05;  % Size of crystal bubbles
phi = rand(nPoints,1);  

Sphi = zeros(nX,1);

for m = 1:nX
    Sphi(m) = ops.S_H(k0, R_b, phi, x(m), 0);
end

plot(x,Sphi,'.')