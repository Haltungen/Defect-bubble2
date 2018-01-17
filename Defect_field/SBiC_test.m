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
nPoints = 10;


R_b = 0.05;  % Size of crystal bubbles
B = shape.Ellipse(R_b, R_b, nPoints);
x = B.points;
phi = rand(nPoints,1);  
Sphi = zeros(nPoints,1);
errorS= zeros(nPoints,1);

S = ops.make_S_K_C(omega, v, v_b, delta, B);
Sphi1 = S*phi;

for m = 1:nPoints
    Sphi(m) = ops.S_BiPeriodic_C(R_b,omega,v,v_b,delta, phi, x(1,m), x(2,m));
    errorS(m) = 100*(Sphi(m)-Sphi1(m))/abs(Sphi1(m));
end 

fprintf('Max error in S: %.5e percentage\n', max(abs(errorS)));