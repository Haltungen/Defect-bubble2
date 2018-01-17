% DemoBubbleResonanceNumericalfunc
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
nPoints = 20;

rho = 5000;
rho_b = 1;
kappa = 5000;
kappa_b = 1;

delta = rho_b/rho;
v = sqrt(rho/kappa);
v_b = sqrt(rho_b/kappa_b);
alpha = [pi, pi];

% NN = 3;     % Order of fourier series
% N1 = 3;     % Order of truncation for lattice sum
% N2 = 5;     % Order of trucation for FFT   
% N3 = 10;    % Number of discretization points in quadrature for integrating G_alpha
R_b = 0.05;  % Size of crystal bubbles
R_d = 0.02;  % Size of defect bubble
B = shape.Ellipse(R_b, R_b, nPoints);
B_d = shape.Ellipse(R_d, R_d, nPoints);

%-------------------------Maximum band-------------------------
initialGuess = 0.259;
z0 = initialGuess;
z1 = initialGuess-initialGuess/100;
z2 = initialGuess-initialGuess/200;    
func_band = @(omega) testfiles.f_Biperiodic(omega, v, v_b, alpha, delta, B); 
char_num_band = tools.MullersMethod(func_band, z2, z1, z0, iterMax, distTol, fTol);

fprintf('Maximum band frequency:                    %.5e  +  i%.5e\n', real(char_num_band), imag(char_num_band));
%Maximum band frequency:                    2.64080e-01  +  i6.65580e-05