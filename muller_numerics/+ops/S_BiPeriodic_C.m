%% S_BiPeriodic_C
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

function out = S_BiPeriodic_C(r,omega,v,v_b,delta, F, x1, x2)
    M = length(x1);
    N = length(F);
    theta = linspace(0,2*pi*(N-1)/N,N).';
    points1 = r*cos(theta);
    points2 = r*sin(theta);
    dsigma = r*2*pi/N;
    Ker = zeros(M,N);
    tol = 10^-6;
    for l = 1:M
        for n = 1:N
            xminusy = sqrt((x1(l)-points1(n))*(x1(l)-points1(n)) + (x2(l)-points2(n))*(x2(l)-points2(n)));
            if xminusy < tol % Singular term
                Ker(l,n) = 1/(2*pi)*(log(dsigma) - 1);
            else
                Ker(l,n) = tools.G(points1(n),points2(n),[x1(l);x2(l)],omega,v,v_b,delta,r);
            end
        end
    end
    out = dsigma*Ker*reshape(F, [], 1);
end