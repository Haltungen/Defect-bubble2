function out  = g(omega, v_w, v_b, delta, radius_b, radius_d, n)
    %fprintf('omega: %.8g    %.8g\n', real(omega), imag(omega));
    
    matrix_A = A_monoPole(omega, v_w, v_b, delta, radius_b, radius_d, n);
% 
%     psi = ones(1, 3);
%     m = matrix_A\psi';
%     
%     s = size(m, 1);
%     denom = abs(dot(m', ones(s, 1)));
%     out = 1/denom;
out = det(matrix_A);
end
