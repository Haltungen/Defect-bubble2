%% f_single
%
% Overview:
%   Evaluates function whose zero correspond to a singular value of the
%   matrix A.
%
% Input:
%   omega:      Frequency
%   v,v_b:      Wave speed in medium and bubble, respectively
%   delta:      Density fraction
%   B:          Bubble representation
%
% Output:
%   out:       The function value

function out  = f_single(omega, v, v_b, delta, B)
    fprintf('Single omega: %.8f    %.8f\n', real(omega), imag(omega));
    nPoints = length(B.points);
    matrix_A = A_single(omega, v, v_b, delta, B);

    psi = ones(1, 2*nPoints);
    m = matrix_A\psi';
    
    s = size(m, 1);
    denom = abs(dot(m', ones(s, 1)));
    out = 1/denom;
end
