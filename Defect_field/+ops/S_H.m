%% S_H
%
% Overview:
%   Evaluates the homogenious Helmholtz single layer potential S_D[F](x).
%
% Input:
%   k:          Wave number
%   r:          Radius of disk
%   F:          Layer density
%   x1,x2:      1:st and 2:nd coordinates for evaluation pts
%
% Output:
%   out:       The single-layer potential (Mx1 vector)

function out = S_H(k, r, F, x1, x2)
    M = length(x1);
    N = length(F);
    theta = linspace(0,2*pi*(N-1)/N,N).';
    points1 = r*cos(theta);
    points2 = r*sin(theta);
    dsigma = r*2*pi/N;
    G = zeros(M,N);
    tol = 10^-6;
    for l = 1:M
        for n = 1:N
            xminusy = sqrt((x1(l)-points1(n))*(x1(l)-points1(n)) + (x2(l)-points2(n))*(x2(l)-points2(n)));
            if xminusy < tol % Singular term
                G(l,n) = 1/(2*pi)*(log(dsigma) - 1);
            else
                G(l,n) = -1i/4*besselh(0,1,(k*xminusy));
            end
        end
    end
    out = dsigma*G*reshape(F, [], 1);
end