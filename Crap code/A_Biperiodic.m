%% A_Biperiodic
%
% Overview:
%   Evaluates the spatial discretization matrix of the boundary integral 
%   operator A.
%
% Input:
%   omega:      Frequency
%   v,v_b:      Wave speed in medium and bubble, respectively
%   delta:      Density fraction
%   B,B_d:      Original and defect bubble, respecitvely
%
% Output:
%   out:        The matrix A

function out = A_Biperiodic(omega, v, v_b, alpha, delta, B)
k = omega*v;
k_b = omega*v_b;
nPoints = length(B.points);

A11 = ops.SingleLayer_H(k_b, B, 'P0', 1);
A12 = ops.SingleLayer_H_Biperiodic(k, B, 'P0', 1, 1, 1, alpha);
A21 = ops.Kstar_H(k_b, B, 'P0', 1);
A22 = ops.Kstar_H_Biperiodic(k, B, 'P0', 1, 1, 1, -alpha);

out = [ A11.Kmat, -1*A12.Kmat; ...
        -0.5*eye(nPoints) + A21.Kmat, -delta*(0.5*eye(nPoints) + A22.Kmat)];
end
