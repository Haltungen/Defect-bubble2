function out = A_freeSpace(omega, v_w, v_b, delta, radius, nPoints)
% A represents the matrix A(omega, delta) in equation (2.3) of Ref 1.
%
% Ref 1: Minnaert resonance for acoustic waves in bubbly media

B = shape.Ellipse(radius, radius, nPoints);

k = omega * v_w;
k_b = omega * v_b;

A11 = ops.SingleLayer_H(k_b, B, 'P0', 1);
A12 = ops.SingleLayer_H(k, B, 'P0', 1);    
A21 = ops.Kstar_H(k_b, B, 'P0', 1);
A22 = ops.Kstar_H(k, B, 'P0', 1);

out = [ A11.Kmat, -1*A12.Kmat ; -1/2.*eye(nPoints) + A21.Kmat, -delta*(1/2.*eye(nPoints) + A22.Kmat) ];

end