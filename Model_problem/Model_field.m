%% main
%
% Overview:
%   We calculate the Minnaert resonance for a single bubble
%   in an infinite extent of liquid by directly solving the
%   boundary integral formulation of the problem for the
%   the characteristic value.
%
% References:
%   Mathematical and Computational Methods in Photonics - Tutorial Notes
%
% Authors:
%   Habib Ammari, Brian Fitzpatrick, Matias Ruiz, Sanghyeon Yu.

clear all
close all

distTol = 5e-4;
fTol = 1e-10;
iterMax = 100;
nPoints = 15;

rho = 5000;
rho_b = 1;
kappa = 5000;
kappa_b = 1;

delta = rho_b/rho;
v = sqrt(rho/kappa);
v_b = sqrt(rho_b/kappa_b);
omega_res = 0.273148681890778 + 1i*0.002988939810069; % Solution from main
eps = real(omega_res)/100;
omega = real(omega_res);
k = omega*v;
k_b = omega*v_b;
r_y = 0;
theta_y = 0;
y = r_y*[cos(theta_y), sin(theta_y)];

% Grid for field
gridMax = 0.03;
gridN = 25;
gridMinX1 = -0;
gridMaxX1 = gridMax;
gridMinX2 = -gridMax;
gridMaxX2 = gridMax;
g1 = linspace(gridMinX1, gridMaxX1, gridN);
g2 = zeros(gridN,1);%linspace(gridMinX2, gridMaxX2, gridN);
%[ g1, g2 ] = meshgrid(g1, g2);
gridPoints = [g1(:) g2(:)]';

gridPointsN = length(gridPoints);

uResponse = zeros(gridPointsN, 1);
uGreen = zeros(gridPointsN, 1);
uTotal = zeros(gridPointsN, 1);

green = @(x1,x2,k) -1i/4*besselh(0,k*sqrt((x1-y(1)).^2+(x2-y(2)).^2));

% NN = 3;     % Order of fourier series
% N1 = 3;     % Order of truncation for lattice sum
% N2 = 5;     % Order of trucation for FFT   
% N4 = 10;    % Number of discretization points in quadrature for integrating G_alphabubble
R_b = 0.05;   % Size of crystal bubbles
R_d = 0.02;  % Size of defect bubble
B = shape.Ellipse(R_b, R_b, nPoints);
B_d = shape.Ellipse(R_d, R_d, nPoints);

RHS = makeRHS(k, k_b, delta, y, B_d);
matrix_A = A(omega, v, v_b, delta, R_b, R_d, nPoints);
Phi = matrix_A\RHS;
phi1 = Phi(1:nPoints);
phi2 = Phi(nPoints+1:2*nPoints);
phi3 = Phi(2*nPoints+1:3*nPoints);

parfor j = 1:gridPointsN
    gridPoint = gridPoints(:, j);    
    % Determine in which region we are in
    if gridPoint(1)^2 + gridPoint(2)^2  <= R_d^2  
        uResponse(j) = ops.S_H(k_b, R_d, phi1, gridPoint(1), gridPoint(2));
        uGreen(j) = green(gridPoint(1), gridPoint(2), k_b);
        uTotal(j) = uResponse(j) + uGreen(j);
    elseif gridPoint(1)^2 + gridPoint(2)^2  <= R_b^2  
        uResponse(j) = ops.S_H(k, R_d, phi2, gridPoint(1), gridPoint(2))+ops.S_H(k, R_b, phi3, gridPoint(1), gridPoint(2));
        uGreen(j) = green(gridPoint(1), gridPoint(2), k);
        uTotal(j) = uResponse(j) + uGreen(j);
    end
end

% Format the fields for plotting
uGreen = reshape(uGreen, [gridN gridN]);
uResponse = reshape(uResponse, [gridN gridN]);
uTotal = reshape(uTotal, [gridN gridN]);

hFig = figure(1);
set(hFig, 'Position', [100 100 1200 900]);
surf(g1, g2, real(uTotal), 'edgecolor', 'none'); xlabel('x1'); ylabel('x2'); title('Greens function for the defect crystal')
axis([gridMinX1, gridMaxX1, gridMinX2, gridMaxX2, min(real(uTotal(:))), max(real(uTotal(:))) ]); rotate3d on;
