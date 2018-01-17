%% G
%
% Overview:
%   Evaluates the Greens function for the bubble crystal, and its normal
%   derivative from outside
%
% Input:
%   x,y:        Points of evaluation (vector, point)
%   omega:      Frequency
%   v,v_b:      Wave speed in medium and bubble, respectively
%   delta:      Density fraction
%   R_b         Radius of bubble
%   NN:         Order of Fourier series
%   N1, N2:     Order of truncation for lattice sum and Fourier coefficient, respectively 
%   N3:         Number of discretization points for integration of S_D^\alpha 
%   N4:         Number of discretization points in quadrature
%
% Output:
%   out:       The Single layer potential

function [out, d_nu_plus] = G(x1,x2,y,omega,v,v_b,delta,R_b,NN,N1,N2,N3,N4)
if nargin < 13 % Default parameters
    NN = 5;
    N1 = 3;
    N2 = 3;     
    N3 = 10;
    N4 = 10;
end
M = length(x1);
out = zeros(M,1);
d_nu_plus = zeros(M,1);
alpha1D = linspace(0,2*pi*(N4-1)/N4,N4);
d2alpha = (1/N4)^2;
parfor j = 1:N4 % Trapezoid quadrature method
    for l = 1:N4
        alpha = [alpha1D(j),alpha1D(l)];
        [g_alpha, g_nu_alpha] = tools.G_alpha(x1,x2,y,omega,v,v_b,alpha,delta,R_b,NN,N1,N2,N3);
        out = out + g_alpha*d2alpha;
        d_nu_plus = d_nu_plus + g_nu_alpha*d2alpha;
    end
end

end