function out  = h(omega, v_w, v_b, delta, radius, nPoints)
    %fprintf('omega: %.8g    %.8g\n', real(omega), imag(omega));
    
    matrix_A = A_freeSpace(omega, v_w, v_b, delta, radius, nPoints);

    psi = ones(1, 2*nPoints);
    m = matrix_A\psi';
    
    s = size(m, 1);
    denom = abs(dot(m', ones(s, 1)));
    out = 1/denom;
end