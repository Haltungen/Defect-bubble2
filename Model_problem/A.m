    function out = A(omega, v_w, v_b, delta, radius_b, radius_d, nPoints)
% A represents the matrix A(omega, delta) in equation (2.3) of Ref 1.
%
% Ref 1: Minnaert resonance for acoustic waves in bubbly media

B = shape.Ellipse(radius_b, radius_b, nPoints);     % Big bubble
B_d = shape.Ellipse(radius_d, radius_d, nPoints);   % Defect bubble

k = omega * v_w;
k_b = omega * v_b;

A11 = ops.SingleLayer_H(k_b, B_d, 'P0', 1);
A12 = ops.SingleLayer_H(k, B_d, 'P0', 1);
A13 = ops.SingleLayer_H(k, B, 'P0', 1, B_d, 'P0', 1);    
A22 = ops.SingleLayer_H(k, B_d, 'P0', 1, B, 'P0', 1);
A23 = ops.SingleLayer_H(k, B, 'P0', 1);
A31 = ops.Kstar_H(k_b, B_d, 'P0', 1);
A32 = ops.Kstar_H(k, B_d, 'P0', 1);
A33 = ops.Kstar_H(k, B, 'P0', 1, B_d, 'P0', 1);


out = [ A11.Kmat, -1.*A12.Kmat, -A13.Kmat ; ...
        zeros(nPoints), 1.*A22.Kmat, A23.Kmat ; ...
        -1/2.*eye(nPoints) + A31.Kmat, -delta.*(1/2.*eye(nPoints) + A32.Kmat), -delta.*A33.Kmat  ];

end

