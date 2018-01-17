%% A
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

function out = A(omega, v, v_b, delta, B, B_d)
k = omega*v;
k_b = omega*v_b;
c = 2/(delta+1);
nPoints = length(B.points);

A11 = ops.SingleLayer_H(k_b, B_d, 'P0', 1);
A12 = ops.SingleLayer_H(k, B_d, 'P0', 1);
A13 = ops.SingleLayer_H(k, B, 'P0', 1, B_d, 'P0', 1);    
A22 = ops.SingleLayer_H(k, B_d, 'P0', 1, B, 'P0', 1);
A23 = ops.SingleLayer_H(k, B, 'P0', 1);
A31 = ops.Kstar_H(k_b, B_d, 'P0', 1);
A32 = ops.Kstar_H(k, B_d, 'P0', 1);
A33 = ops.Kstar_H(k, B, 'P0', 1, B_d, 'P0', 1);
A42 = ops.Kstar_H(k, B_d, 'P0', 1, B, 'P0', 1);
A43 = ops.Kstar_H(k, B, 'P0', 1);
[A24, A44] = ops.make_S_K_C(omega, v, v_b, delta, B);

out = [ A11.Kmat, -1*A12.Kmat, -1*A13.Kmat, zeros(nPoints); ...
        zeros(nPoints), A22.Kmat, A23.Kmat, -1*A24; ...
        -1/2.*eye(nPoints) + A31.Kmat, -delta*(1/2.*eye(nPoints) + A32.Kmat), -delta*A33.Kmat, zeros(nPoints); ...
        zeros(nPoints), A42.Kmat, -1/2.*eye(nPoints) + A43.Kmat, -1*(c*1/2.*eye(nPoints)+ A44)];

end

