%% KstarBiPeriodic
%
% Overview:
%   Evaluates the quasi-biperiodic Neumann-Poincare operator K_D^*[F](x) of
%   a disk D.
%
% Input:
%   k:          Wave number
%   alpha:      Quasi-periodicity
%   r:          Radius of disk
%   F:          Layer density
%   x1,x2:      1:st and 2:nd coordinates for evaluation pts
%
% Output:
%   out:       The Neumann-Poincare operator (Mx1 vector)

function out = KstarBiPeriodic(k, alpha, r, F, x1, x2) 
    M = length(x1);
    N = length(F);
    theta = linspace(0,2*pi*(N-1)/N,N);
    points1 = r*cos(theta);
    points2 = r*sin(theta);
    dsigma = r*2*pi/N;
    G = zeros(M,N);
    tol = (10^-6)^2;
    for l = 1:M
        for n = 1:N
            if (points1(n)-x1(l))^2+(points2(n)-x2(l))^2 < tol % Singular term
                G(l,n) = 1/(4*pi*r);
            else
                theta_x = atan2(x2(l),x1(l));
                G(l,n) = dot([cos(theta_x);sin(theta_x)],ops.GradGBiPeriodic(k, x1(l)-points1(n), x2(l)-points2(n), 1, 1, -alpha));
            end
        end
    end
    out = dsigma*G*reshape(F, [], 1);
end