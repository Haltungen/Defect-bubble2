rho = 5000;
rho_b = 1;
kappa = 5000;
kappa_b = 1;

delta = rho_b/rho;
v = sqrt(rho/kappa);
v_b = sqrt(rho_b/kappa_b)*0.5;
omega = 1;
k = omega*v;
k_b = omega*v_b;
r_y = 0;
theta_y = 0;
y = r_y*[cos(theta_y), sin(theta_y)];
c = 2/(delta+1);

NN = 3;
N1 = 1;
N2 = 1;
N3 = 5;
N4 = 5;
nPoints = 10000;

R_b = 0.05;   % Size of crystal bubbles
B = shape.Ellipse(R_b, R_b, nPoints);
alpha = [pi/2,pi];

sigma_B = B.sigma;
pts1 = B.points(1,:);
pts2 = B.points(2,:);
i=1;
[g1,g1_nu] = tools.G_alpha(pts1(i),pts2(i),[pts1(i+1),pts2(i+1)],omega,v,v_b,alpha,delta,R_b,NN,N1,N2,N3);
S = g1*sigma_B(i);
Sref = c*sigma_B(i)/(2*pi)*(log(sigma_B(i)) - 1);
K = g1_nu*sigma_B(i);
Kref = sigma_B(i)*c/(4*pi*R_b);
Green0 = tools.G_alpha(0,eps,[0,0],omega,v,v_b,alpha,delta,R_b,NN,N1,N2,N3);
Green0ref = c/(2*pi)*log(eps);
fprintf('Relative error in diagonal of S: %0.5f\n', (S-Sref)/Sref)
fprintf('Relative error in diagonal of K: %0.5f\n', (K-Kref)/Kref);
fprintf('Relative error in approx. of Greens function singularity: %0.5f\n', (Green0-Green0ref)/Green0ref);
