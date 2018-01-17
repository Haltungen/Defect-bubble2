    classdef Kstar_H_Biperiodic < ops.Operators
    
    methods
        function obj = Kstar_H_Biperiodic(k, D, type, step, d1, d2, alpha)
            obj = obj@ops.Operators(D, type, step, D, type, step);
            
            obj.Kmat = ops.Kstar_H_Biperiodic.make_kernel_matrix(k, D.points, D.tvec, D.normal, D.avec, D.sigma, d1, d2, alpha);
        end
        
    end
    
    methods(Static)
        function K = make_kernel_matrix(k, D, tvec, normal, avec, sigma, d1, d2, alpha)
            N = size(D,2);
            K = zeros(N, N);
             
            for i = 1:N
                grad = tools.GradGBiPeriodic(k, D(1,i)-D(1,:), D(2,i)-D(2,:), d1, d2, -alpha);
                K(i, :) = (normal(:,i)'*grad).*sigma;
                
                tvec_norm_square = tvec(1,:).^2 + tvec(2,:).^2;
                
                if i == 1
                    K(i,i) = 1/2 * (K(i, N) +  K(i, 2));
                elseif i == N
                    K(i,i) = 1/2 * (K(i, i-1) +  K(i, 1));
                else
                   K(i,i) = -1/4/pi*avec(:,i)'*normal(:,i)/tvec_norm_square(i)*sigma(i);
%                    K(i,i) = 1/2 * (K(i, i-1) +  K(i, i+1));
                end                
                
            end
            %%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        function [ val ] = eval(D, F)
            % Evaluate K_D^*[F]
            % INPUTS:
            % D: shape, a C2boundary object
            % F: function defined on D
            
            Ks = ops.Kstar.make_kernel_matrix(D.points, D.tvec, D.normal, D.avec, ...
                D.sigma);
            val = Ks*F(:);
        end
        
    end
end

