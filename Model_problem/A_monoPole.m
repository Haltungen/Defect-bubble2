function out = A_monoPole(omega, v_w, v_b, delta, r_b, r_d, n)

k = omega*v_w;
k_b = omega*v_b;
out = [besselj(n,k_b*r_d), -besselj(n,k*r_d), -besselh(n,1,k*r_d); ...
       0.5*k_b*(besselj(n-1,k_b*r_d)-besselj(n+1,k_b*r_d)), -delta*0.5*k*(besselj(n-1,k*r_d)-besselj(n+1,k*r_d)), -delta*k*0.5*(besselh(n-1,1,k*r_d)-besselh(n+1,1,k*r_d)); ...
       0, besselj(n,k*r_b), besselh(n,1,k*r_b)];
end
