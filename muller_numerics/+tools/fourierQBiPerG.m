%% fourierQBiPerG
%
% Overview:
%   Returns the Fourier coefficients of the quasi-biperiodic function \Gamma_\alpha^k(x,y) as as a function of x on a circle of radius r.
%
% Input:
%   r:          Radius of the circle containing x
%   y:          y
%   alp:        \alpha
%   k:          k
%   NN:         Order of Fourier series
%   N1:         Order of the sum when computing the lattice sums
%   N2:         Order of the sum for the Fourier coefficients
%
% Output:
%   out:       The Fourier coefficients.
%
% References:
%   Mathematical and Computational Methods in Photonics - Tutorial Notes
%
% Authors:
%   Habib Ammari, Brian Fitzpatrick, Sanghyeon Yu, Erik Orvehed Hiltunen

function [out, out_nu] = fourierQBiPerG(r, y, alp, k, NN, N1, N2) 

%%% store lattice sum data %%%
L = 1;
data_Qn_posi = zeros(NN+N2+1,1);
data_Qn_nega = zeros(NN+N2,1);
for j=1:N2+NN+1
    data_Qn_posi(j)=tools.lattice_Sn(j-1,k,alp,L,N1);
end
for j=1:N2+NN   
    data_Qn_nega(j)=tools.lattice_Sn(j,k,[-alp(1), alp(2)],L,N1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r_2 = sqrt(y(1)*y(1)+y(2)*y(2));
theta_2 = atan2(y(2),y(1));
out = zeros(1,2*NN+1);
out_nu = zeros(1,2*NN+1);
for n = -NN:NN
    sgn = (-1)^mod(n,2);
    term = 0;
    for l = -N2:N2
        if -n-l >= 0
            Qn = data_Qn_posi(-(n+l)+1);  
        else
            Qn = data_Qn_nega(n+l);
        end
        term = term + Qn*besselj(l,k*r_2)*exp(1i*l*theta_2);
    end
    out(NN+1+n) = term*sgn*besselj(n,k*r);
    out_nu(NN+1+n) = term*sgn*0.5*k*(besselj(n-1,k*r) - besselj(n+1,k*r));
    if abs(r-r_2) < 10^-6
        out(NN+1+n) = out(NN+1+n) + sgn*besselh(n,1,k*r)*besselj(-n,k*r_2)*exp(-1i*n*theta_2);
        out_nu_plus = 0.5*k*sgn*besselh(-n,k*r_2)*exp(-1i*n*theta_2)*(besselj(n-1,k*r) - besselj(n+1,k*r));
        out_nu_minus = 0.5*k*sgn*(besselh(n-1,k*r) - besselh(n+1,k*r))*besselj(-n,k*r_2)*exp(-1i*n*theta_2); 
        out_nu(NN+1+n) = out_nu(NN+1+n) + 0.5*(out_nu_minus+out_nu_plus);     
    elseif r_2 > r
        out(NN+1+n) = out(NN+1+n) + sgn*besselh(-n,1,k*r_2)*exp(-1i*n*theta_2)*besselj(n,k*r);
        out_nu(NN+1+n) = out_nu(NN+1+n) + 0.5*k*sgn*besselh(-n,k*r_2)*exp(-1i*n*theta_2)*(besselj(n-1,k*r) - besselj(n+1,k*r));
    else
        out(NN+1+n) = out(NN+1+n) + sgn*besselh(n,1,k*r)*besselj(-n,k*r_2)*exp(-1i*n*theta_2);
        out_nu(NN+1+n) = out_nu(NN+1+n) + 0.5*k*sgn*(besselh(n-1,k*r) - besselh(n+1,k*r))*besselj(-n,k*r_2)*exp(-1i*n*theta_2);
    end
end
out = -1i/4.*out.';
out_nu = -1i/4.*out_nu.';
end