function [ G GSpec GSpat ] = GBiPeriodic(k, x1, x2, d1, d2, alpha)

E = pi/d1/d2;

%%%%%%%%%%%%%%%%%%%%%%%%%
% GSpec
%%%%%%%%%%%%%%%%%%%%%%%%%
pqMax = 1;

GSpec = zeros(size(x1));
parfor p=-pqMax:pqMax
    for q=-pqMax:pqMax
        kxp = -alpha(1) + 2*pi*p/d1;
        kyq = -alpha(2) + 2*pi*q/d2;
        
        kpqSqr = kxp.^2 + kyq.^2;
        gammapqSqr = kpqSqr - k^2;
        
        GSpec = GSpec + exp(-gammapqSqr/4/E).*exp(-1i*(kxp*x1 + kyq*x2))./gammapqSqr;
    end
end

GSpec = -1/d1/d2 * GSpec;

%%%%%%%%%%%%%%%%%%%%%%%%%
% GSpat
%%%%%%%%%%%%%%%%%%%%%%%%%
mnMax = 1;
Q = 3;

GSpat = 0;
for m=-mnMax:mnMax
    for n=-mnMax:mnMax
        rhomn = [ m*d1; n*d2 ];
        rho_rhomnSqr = (x1-rhomn(1)).^2 + (x2-rhomn(2)).^2;
        
        term1 = 0;
        Eqp1 = ops.E1(rho_rhomnSqr*E);
        for q=0:Q
            term1 = term1 + (k/2/sqrt(E)).^(2*q)./factorial(q).*Eqp1;
            Eqp1 = (exp(-rho_rhomnSqr*E) - rho_rhomnSqr*E.*Eqp1)./(q+1);
        end
        
        GSpat = GSpat - 1/4/pi.*exp(1i*(alpha(1)*rhomn(1) + alpha(2)*rhomn(2))).*term1;
    end
end

G = GSpat + GSpec;

end