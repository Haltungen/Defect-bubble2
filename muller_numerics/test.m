nPoints = 10;
R_b = 0.05;   % Size of crystal bubbles
R_d = 0.045;  % Size of defect bubble
k = 1;
B = shape.Ellipse(R_b, R_b, nPoints);
B_d = shape.Ellipse(R_d, R_d, nPoints);

A21 = ops.Kstar_H(k, B, 'P0', 1);
A21.Kmat