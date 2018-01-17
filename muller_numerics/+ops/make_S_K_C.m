function [outS, outK] = make_S_K_C(omega, v, v_b, delta, D)
%% make_S_K_C
%
% Overview:
%   Returns the kernel matrices of S_D^\# and K_D^\#. These computations
%   are combined in one function instead of two, to avoid recomputation
%   of the same values. 
%
% Input:
%   r:          Radius of the circle containing x
%   y:          y
%   alp:        \alpha
%   k:          k
%   NN:         Limits for the Fourier coefficients, between -NN and NN
%   N:          Order of the sum when computing the lattice sums
%   M:          Order of the sum for the Fourier coefficients
%
% Output:
%   outS:       Kernel matrix of S_D^\#
%   outK:       Kernel matrix of K_D^\#
%
% References:
%   Mathematical and Computational Methods in Photonics - Tutorial Notes
%
% Authors:
%   Habib Ammari, Brian Fitzpatrick, Sanghyeon Yu, Erik Orvehed Hiltunen

pts = D.points;
sigma_D = D.sigma;
R_b = D.axis_a;
green = @(x1,x2,y) tools.G(x1,x2,y,omega,v,v_b,delta,R_b);

N = size(pts,2);
outS = zeros(N);
outK = zeros(N);

for i=1:N
    [g1, g_nu1] = green(pts(1,1:(i-1)),pts(2,1:(i-1)),[pts(1,i),pts(2,i)]);
    [g2, g_nu2] = green(pts(1,(i+1):end),pts(2,(i+1):end),[pts(1,i),pts(2,i)]);
    outS(1:(i-1),i) = sigma_D(i)*g1;
    outS((i+1):end,i) = sigma_D(i)*g2;
    outK(1:(i-1),i) = sigma_D(i)*g_nu1;
    outK((i+1):end,i) = sigma_D(i)*g_nu2;

    % Diagonal terms.
    c = 2/(delta+1);
    outS(i,i) = c/(2*pi)*sigma_D(i)*(log(sigma_D(i)) - 1);
    outK(i,i) = sigma_D(i)*c/(4*pi*R_b); % page 22 in Tutorial notes, for D=circle
    %outS(i,i) = 1/2 * (outS(i, i-1) +  outS(i, i+1));
    %outK(i,i) = 1/2 * (outK(i, i-1) +  outK(i, i+1));
end
end