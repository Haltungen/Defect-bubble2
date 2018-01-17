%% A_single
%
% Overview:
%   Evaluates the spatial discretization matrix of the boundary integral 
%   operator A for the case of a single bubble in free space.
%
% Input:
%   omega:      Frequency
%   v,v_b:      Wave speed in medium and bubble, respectively
%   delta:      Density fraction
%   B:          Bubble representation
%
% Output:
%   out:        The matrix A

function out = A_single(omega, v, v_b, delta, B)
nPoints = length(B.points);
k = omega * v;
k_b = omega * v_b;

A11 = ops.SingleLayer_H(k_b, B, 'P0', 1);
A12 = ops.SingleLayer_H(k, B, 'P0', 1);    
A21 = ops.Kstar_H(k_b, B, 'P0', 1);
A22 = ops.Kstar_H(k, B, 'P0', 1);

out = [ A11.Kmat, -1*A12.Kmat ; -1/2.*eye(nPoints) + A21.Kmat, -delta*(1/2.*eye(nPoints) + A22.Kmat) ];

end

