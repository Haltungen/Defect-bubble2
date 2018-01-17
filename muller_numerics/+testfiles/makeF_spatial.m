%% makeF
%
% Overview:
%   Evaluates the RHS F in the boundary integral equation. Also returns the
%   Fourier coefficients of the Greens function \Gamma.
%
% Input:
%   r:          Radius of bubble
%   y:          Points of evaluation (vector, point)
%   omega:      Frequency
%   v,v_b:      Wave speed in medium and bubble, respectively
%   alpha:      Quasi-periodicty
%   delta:      Density fraction
%   NN:         Order of Fourier series
%   N:          Order of the sum when computing the lattice sums
%   M:          Order of the sum for the Fourier coefficients
%
% Output:
%   out:        The Single layer potential
%   Gamma_k:    Fourier coefficients of \Gamma_k (and kb)
%   Gamma_k_nu: Fourier coefficients of \partial_\nu \Gamma_k (and kb)

function out = makeF_spatial(B, y, omega, v, v_b, alpha, delta)
k = omega*v;
k_b = omega*v_b;
x1 = B.points(1,:);
x2 = B.points(2,:);
theta = atan2(x2,x1);

Gamma_k = tools.GBiPeriodic(k, x1-y(1), x2-y(2), 1, 1, alpha);
Gamma_kb = tools.GBiPeriodic(k_b, x1-y(1), x2-y(2), 1, 1, alpha);
Gamma_k_nu = dot([cos(theta);sin(theta)],tools.GradGBiPeriodic(k, x1-y(1), x2-y(2), 1, 1, alpha));
Gamma_kb_nu = dot([cos(theta);sin(theta)],tools.GradGBiPeriodic(k_b, x1-y(1), x2-y(2), 1, 1, alpha));
F1 = Gamma_k - Gamma_kb;
F2 = -delta*Gamma_k_nu - Gamma_kb_nu;
out = [F1.'; F2.'];

end