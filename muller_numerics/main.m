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

distTol = 5e-3;
fTol = 5e-3;
iterMax = 100;
nPoints = 20;

rho = 5000;
rho_b = 1;
kappa = 5000;
kappa_b = 1;

delta = rho_b/rho;
v = sqrt(rho/kappa);
v_b = sqrt(rho_b/kappa_b);


% NN = 3;     % Order of fourier series
% N1 = 3;     % Order of truncation for lattice sum
% N2 = 5;     % Order of trucation for FFT   
% N4 = 10;    % Number of discretization points in quadrature for integrating G_alpha
R_b = 0.05;   % Size of crystal bubbles
R_d = 0.045;  % Size of defect bubble
B = shape.Ellipse(R_b, R_b, nPoints);
B_d = shape.Ellipse(R_d, R_d, nPoints);

%-------------------------Single bubble-------------------------
initialGuess = 0.0002;
z0 = initialGuess;
z1 = initialGuess-initialGuess/100;
z2 = initialGuess-initialGuess/200;    
func_single = @(omega) f_single(omega, v, v_b, delta, B);
char_num_single = tools.MullersMethod(func_single, z2, z1, z0, iterMax, distTol, fTol);


%-------------------------Maximum band-------------------------
initialGuess = 0.259;
z0 = initialGuess;
z1 = initialGuess-initialGuess/100;
z2 = initialGuess-initialGuess/200;    
func_band = @(omega) f_band(omega, v, v_b, delta, R_b);
char_num_band = tools.MullersMethod(func_band, z2, z1, z0, iterMax, distTol, fTol);

fprintf('Frequency of defect bubble in free space:  %.5e  +  i%.5e\n', real(char_num_single), imag(char_num_single));
fprintf('Maximum band frequency:                    %.5e  +  i%.5e\n', real(char_num_band), imag(char_num_band));

%-------------------------Defect bubble-------------------------
eps = (abs(real(char_num_band)) + imag(char_num_band))/10; % Start inside the band-gap
initialGuess = abs(real(char_num_band)) + imag(char_num_band) + eps;
initialGuess = 0.27314586 -1i*0.00298593;
z0 = initialGuess;
z1 = initialGuess-initialGuess/100;
z2 = initialGuess-initialGuess/200;    
func = @(omega) f(omega, v, v_b, delta, B, B_d);
char_num = tools.MullersMethod(func, z2, z1, z0, iterMax, distTol, fTol);
fprintf('Frequency of defect bubble in free space:  %.5e  +  i%.5e\n', real(char_num_single), imag(char_num_single));
fprintf('Maximum band frequency:                    %.5e  +  i%.5e\n', real(char_num_band), imag(char_num_band));
fprintf('Frequency close to the defect bubble:      %.5e  +  i%.5e\n', real(char_num), imag(char_num));

%% nPoints = 20, R_d = 0.02, Defect omega: 0.63309181    -0.00086662

