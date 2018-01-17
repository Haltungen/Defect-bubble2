function val = wavenb_acs(freq, rho, lambda)
% Wavenumber
% rho, lambda: density and bulk modulus
    val = freq*sqrt(rho/lambda);
end