clear all
close all

rho = 5000;
rho_b = 1;
kappa = 5000;
kappa_b = 5000;
omega = 0.27;

delta = rho_b/rho;
v = sqrt(rho/kappa);
v_b = sqrt(rho_b/kappa_b);

k0 = omega * v;
kb = omega * v_b;

NN = 2;     % Order of fourier series
N1 = 1;     % Order of truncation for lattice sum
N2 = 1;     % Order of trucation inside Fourier coefficients  
N3 = 5;     % Number of discretization points for integration of S_D^\alpha 
R_b = 0.05;
alpha = [pi/2,pi];
nPoints = 10;

xMax = 0.0501;
xMin = 0.0499;
nX = 30;
hX = (xMax-xMin)/(nX-1);
x1 = xMin:hX:xMax;
r_2 = R_b;
theta_2 = pi/2+0.00003;
y = r_2*[cos(theta_2),sin(theta_2)];
B = shape.Ellipse(R_b, R_b, nPoints);
Green = tools.G_alpha(x1,zeros(1,nX),y,omega,v,v_b,alpha,delta,R_b,NN,N1,N2,N3);
Green_spatial = testfiles.G_alpha_spatial(x1,zeros(1,nX),y,omega,v,v_b,alpha,delta,B);
figure
plot(x1,real(Green),'.',x1,real(Green_spatial),'.');
legend('Spectral','Spatial')
figure
plot(x1,imag(Green),'.',x1,imag(Green_spatial),'.');
legend('Spectral','Spatial')
fprintf('Maximal error: %.5e\n',max(abs(Green-Green_spatial)));