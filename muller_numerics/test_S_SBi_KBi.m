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
xMax = 0.5;
nX = 3^2;
hX = 2*xMax/(nX-1);
[x1, x2] = meshgrid(-xMax:hX:xMax);
alpha = [pi/2, pi/8];

R_b = 0.05;  % Size of crystal bubbles
R_d = 0.02;  % Size of defect bubble
B = shape.Ellipse(R_b, R_b, nPoints);
B_d = shape.Ellipse(R_d, R_d, nPoints);
x = B.points;
phi = rand(nPoints,1);  
Sbiphi = zeros(nPoints,1);
Kphi = zeros(nPoints,1);
Sphi = zeros(nPoints,1);
errorSbi= zeros(nPoints,1);
errorS= zeros(nPoints,1);
errorK = zeros(nPoints,1);

Sbi = ops.SingleLayer_H_Biperiodic(k0, B, 'P0', 1, 1, 1, alpha);
K = ops.Kstar_H_Biperiodic(k0, B, 'P0', 1, 1, 1, alpha);
S = ops.SingleLayer_H(k0, B, 'P0', 1);
Sbiphi1 = Sbi.Kmat*phi;
Kphi1 = K.Kmat*phi;
Sphi1 = S.Kmat*phi;

for m = 1:nPoints
    Sbiphi(m) = ops.SBiPeriodic(k0, alpha, R_b, phi, x(1,m), x(2,m));
    Kphi(m) = ops.KstarBiPeriodic(k0, alpha, R_b, phi, x(1,m), x(2,m));
    Sphi(m) = ops.S_H(k0, R_b, phi, x(1,m), x(2,m));
    errorS(m) = (Sphi(m)-Sphi1(m))/abs(Sphi1(m));
    errorSbi(m) = (Sbiphi(m)-Sbiphi1(m))/abs(Sbiphi1(m));
    errorK(m) = (Kphi(m)-Kphi1(m))/abs(Kphi1(m));
end

fprintf('Max error in S: %.5e\n', max(abs(errorS)));
fprintf('Max error in SBi: %.5e\n', max(abs(errorSbi)));
fprintf('Max error in KBi: %.5e\n', max(abs(errorK)));