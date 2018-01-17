function [ GradG GradGSpec GradGSpat ] = GradGBiPeriodic(k, x1, x2, d1, d2, alpha)

E = pi/d1/d2;

%%%%%%%%%%%%%%%%%%%%%%%%%
% GSpec
%%%%%%%%%%%%%%%%%%%%%%%%%
pqMax = 1;

GradGSpec1 = zeros(size(x1));
GradGSpec2 = zeros(size(x1));
parfor p=-pqMax:pqMax
    for q=-pqMax:pqMax
        kxp = -alpha(1) + 2*pi*p/d1;
        kyq = -alpha(2) + 2*pi*q/d2;
        
        kpqSqr = kxp.^2 + kyq.^2;
        gammapqSqr = kpqSqr - k^2;
        
        GradGSpec1 = GradGSpec1 + 1i*kxp*exp(-gammapqSqr/4/E).*exp(-1i*(kxp*x1 + kyq*x2))./gammapqSqr;
        GradGSpec2 = GradGSpec2 + 1i*kyq*exp(-gammapqSqr/4/E).*exp(-1i*(kxp*x1 + kyq*x2))./gammapqSqr;
    end
end

GradGSpec1 = 1/d1/d2 * GradGSpec1;
GradGSpec2 = 1/d1/d2 * GradGSpec2;

GradGSpec = [ GradGSpec1 ; GradGSpec2 ];

%%%%%%%%%%%%%%%%%%%%%%%%%
% GSpat
%%%%%%%%%%%%%%%%%%%%%%%%%
mnMax = 1;
Q = 3;

GradGSpat1 = 0;
GradGSpat2 = 0;
parfor m=-mnMax:mnMax
    for n=-mnMax:mnMax
        rhomn = [ m*d1; n*d2 ];
        rho_rhomnSqr = (x1-rhomn(1)).^2 + (x2-rhomn(2)).^2;
        
        term1 = E/2/pi*(x1-rhomn(1))*exp(1i*(alpha(1)*rhomn(1) + alpha(2)*rhomn(2)));
        term2 = E/2/pi*(x2-rhomn(2))*exp(1i*(alpha(1)*rhomn(1) + alpha(2)*rhomn(2)));
           
        term3 = exp(-rho_rhomnSqr*E)./(rho_rhomnSqr*E);
        Eq = ops.E1(rho_rhomnSqr*E);
        for q=1:Q
            term3 = term3 + (k/2/sqrt(E)).^(2*q)./factorial(q).*Eq;
            Eq = (exp(-rho_rhomnSqr*E) - rho_rhomnSqr*E.*Eq)./q;%Eqp1 = (exp(-rho_rhomnSqr*E) - rho_rhomnSqr*E.*Eqp1)./(q+1);
        end
        
        GradGSpat1 = GradGSpat1 + term1.*term3;
        GradGSpat2 = GradGSpat2 + term2.*term3;
    end
end

GradGSpat = [ GradGSpat1 ; GradGSpat2 ];

GradG = GradGSpat + GradGSpec;

end

%%
% function [ GradG GradGSpec GradGSpat ] = GradGBiperiodic(k, x1, x2, d1, d2, alpha)
% 
% E = pi/d1/d2;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % GSpec
% %%%%%%%%%%%%%%%%%%%%%%%%%
% pqMax = 1;
% 
% GradGSpec1 = zeros(size(x1));
% GradGSpec2 = zeros(size(x1));
% parfor p=-pqMax:pqMax
%     for q=-pqMax:pqMax
%         kxp = -alpha(1) + 2*pi*p/d1;
%         kyq = -alpha(2) + 2*pi*q/d2;
%         
%         kpqSqr = kxp.^2 + kyq.^2;
%         gammapqSqr = kpqSqr - k^2;
%         
%         GradGSpec1 = GradGSpec1 + 1i*kxp*exp(-gammapqSqr/4/E).*exp(-1i*(kxp*x1 + kyq*x2))./gammapqSqr;
%         GradGSpec2 = GradGSpec2 + 1i*kyq*exp(-gammapqSqr/4/E).*exp(-1i*(kxp*x1 + kyq*x2))./gammapqSqr;
%     end
% end
% 
% GradGSpec1 = 1/d1/d2 * GradGSpec1;
% GradGSpec2 = 1/d1/d2 * GradGSpec2;
% 
% GradGSpec = [ GradGSpec1 ; GradGSpec2 ];
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%
% % GSpat
% %%%%%%%%%%%%%%%%%%%%%%%%%
% mnMax = 1;
% Q = 3;
% 
% GradGSpat1 = 0;
% GradGSpat2 = 0;
% parfor m=-mnMax:mnMax
%     for n=-mnMax:mnMax
%         rhomn = [ m*d1; n*d2 ];
%         rho_rhomnSqr = (x1-rhomn(1)).^2 + (x2-rhomn(2)).^2;
%         
%         term1 = E/2/pi*(x1-rhomn(1))*exp(1i*(alpha(1)*rhomn(1) + alpha(2)*rhomn(2)));
%         term2 = E/2/pi*(x2-rhomn(2))*exp(1i*(alpha(1)*rhomn(1) + alpha(2)*rhomn(2)));
%            
%         term3 = exp(-rho_rhomnSqr*E)./(rho_rhomnSqr*E);
%         Eq = ops.E1(rho_rhomnSqr*E);
%         for q=1:Q
%             term3 = term3 + (k/2/sqrt(E)).^(2*q)./factorial(q).*Eq;
%             Eq = (exp(-rho_rhomnSqr*E) - rho_rhomnSqr*E.*Eq)./q;
%         end
%         
%         GradGSpat1 = GradGSpat1 + term1.*term3;
%         GradGSpat2 = GradGSpat2 + term2.*term3;
%     end
% end
% 
% GradGSpat = [ GradGSpat1 ; GradGSpat2 ];
% 
% GradG = GradGSpat + GradGSpec;
% 
% end