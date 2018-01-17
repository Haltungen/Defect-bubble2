%% f
%
% Overview:
%   Evaluates function whose zero correspond to a singular value of the
%   matrix A.
%
% Input:
%   omega:      Frequency
%   v,v_b:      Wave speed in medium and bubble, respectively
%   delta:      Density fraction
%   B,B_d:      Original and defect bubble, respecitvely
%
% Output:
%   out:       The function value

function out  = f(omega, v, v_b, delta, B, B_d)
    fprintf('Defect omega: %.8f    %.8f\n', real(omega), imag(omega));
    
    matrix_A = A(omega, v, v_b, delta, B, B_d);
    
    psi = ones(1, 4*length(B.points));
    m = matrix_A\psi';
    
    s = size(m, 1);
    denom = abs(dot(m', ones(s, 1)));
    out = 1/denom;
end
