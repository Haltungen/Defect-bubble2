function out  = f(omega, v_w, v_b, delta, radius_b, radius_d, nPoints)
    %fprintf('omega: %.8g    %.8g\n', real(omega), imag(omega));
    
    matrix_A = A(omega, v_w, v_b, delta, radius_b, radius_d, nPoints);

    psi = ones(1, 3*nPoints);
    m = matrix_A\psi';
    
    s = size(m, 1);
    denom = abs(dot(m', ones(s, 1)));
    out = 1/denom;
end
