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


function [out, Gamma_k, Gamma_kb, Gamma_k_nu, Gamma_kb_nu] = makeF(r, y, omega, v, v_b, alpha, delta, NN, N, M)
k = omega*v;
k_b = omega*v_b;

[Gamma_k,Gamma_k_nu] = tools.fourierQBiPerG(r, y, alpha, k, NN, N, M);
[Gamma_kb, Gamma_kb_nu] = tools.fourierQBiPerG(r, y, alpha, k_b, NN, N, M);
F1 = Gamma_k - Gamma_kb;
F2 = -delta*Gamma_k_nu - Gamma_kb_nu;
out = [F1; F2];

end