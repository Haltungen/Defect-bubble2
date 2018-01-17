%% G_alpha_spatial
%
% Overview:
%   Evaluates the quasi-biperiodic Greens function for the bubble crystal,
%   and its normal derivative from outside
%
% Input:
%   x,y:        Points of evaluation (vector, point)
%   omega:      Frequency
%   v,v_b:      Wave speed in medium and bubble, respectively
%   alpha:      Quasi-periodicity
%   delta:      Density fraction
%   R_b         Radius of bubble
%   NN:         Order of Fourier series
%   N1, N2:     Order of truncation for lattice sum and Fourier series, respectively 
%   N3:         Number of discretazation points for integration of S_D^\alpha 
%
% Output:
%   out:        The Green's function (Mx1 vector)
%   d_nu_p:     The radial derivative for x->\partial D from
%               outside.
                

function [out, d_nu_plus] = G_alpha_spatial(x1,x2,y,omega,v,v_b,alpha,delta,B)
k = omega*v;
k_b = omega*v_b;
L = 1;

A = A_Biperiodic(omega, v, v_b, alpha, delta, B);
F = makeF_spatial(B, y, omega, v, v_b, alpha, delta);
R_b = B.axis_a;

Phi = A\F;
phi_b = Phi(1:end/2);
phi = Phi(end/2+1:end);

M = length(x1);
out = zeros(M,1);
d_nu_plus = zeros(M,1);
for j = 1:M
    if sqrt(x1(j)*x1(j)+x2(j)*x2(j)) < R_b 
        out(j) = ops.S_H(k_b, R_b, phi_b, x1(j), x2(j))+ops.GBiPeriodic(k_b, x1(j)-y(1), x2(j)-y(2), L, L, alpha);
    else
        out(j) = ops.SBiPeriodic(k, alpha, R_b, phi, x1(j), x2(j))+ops.GBiPeriodic(k, x1(j)-y(1), x2(j)-y(2), L, L, alpha);
    end
    d_nu_plus(j) = NaN;%sum((Gamma_k_nu+0.5*phi_hat).*exp(1i.*n.*theta_x)) + ops.KstarBiPeriodic(k, alpha, R_b, phi, x1(j), x2(j)); % Assumes x-> \partial D +
end