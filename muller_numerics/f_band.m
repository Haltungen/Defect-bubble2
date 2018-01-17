%% f_band
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

function out  = f_band(omega, v, v_b, delta, R)
    fprintf('Band omega: %.8f    %.8f\n', real(omega), imag(omega));
    NN = 5;
    N = 3;
    matrix_A = MakeA_Bubble(omega,v,v_b,[pi,pi],delta,1,R,NN,N);
    psi = ones(1, 2*(2*NN+1));
    m = matrix_A\psi';
    
    s = size(m, 1);
    denom = abs(dot(m', ones(s, 1)));
    out = 1/denom;
end
