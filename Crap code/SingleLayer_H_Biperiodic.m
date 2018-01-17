classdef SingleLayer_H_Biperiodic < ops.Operators
    
    methods
        function obj = SingleLayer_H_Biperiodic(k, D1, type1, step1, d1, d2, alpha)
            D2 = D1;
            type2 = type1;
            step2 = step1;
            
            obj = obj@ops.Operators(D1, type1, step1, D2, type2, step2);
            obj.Kmat = ops.SingleLayer_H_Biperiodic.make_kernel_matrix(k, D1.points, D1.sigma, D1.tvec, d1, d2, alpha);
        end
        
    end
    
    methods(Static)
        function K = make_kernel_matrix(k, D, sigma_D, tvec, d1, d2, alpha)            
            N = size(D,2);
            K = zeros(N) ;
                       
            for i=1:N
                Gi = ops.GBiperiodic(k, D(1,i)-D(1,:), D(2,i)-D(2,:), d1, d2, alpha);
                K(i,:) = sigma_D.*Gi;

                 if i == 1
                     K(i,i) = 1/2 * (K(i, N) +  K(i, 2));
                 elseif i == N
                     K(i,i) = 1/2 * (K(i, i-1) +  K(i, 1));
                 else
%                     E = pi/d1/d2;
%                      K(i,i) = 1/2/pi*sigma_D(i)*(log(sqrt(E)*sigma_D(i)) - 1);
                    K(i,i) = 1/2 * (K(i, i-1) +  K(i, i+1));
                 end
            end
        end
        
        function [ val ] = eval(D, F, X)
        end
    end
    
end

