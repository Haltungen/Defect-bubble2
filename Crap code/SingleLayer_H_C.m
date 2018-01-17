%% Class for Single Layer potential operator (Helmholtz equation)
% Single layer potential for bubble crystal medium
%
classdef SingleLayer_H_C < ops.Operators
	% Class for single layer potential operator for Helmholtz equation
	
	methods
		function obj = SingleLayer_H_C(omega, v, v_b, delta, D1, type1, step1)
            
			obj = obj@ops.Operators(D1, type1, step1);
            obj.Kmat = ops.SingleLayer_H.make_kernel_matrix(omega, v, v_b, delta, D1.points, D1.sigma);
		end
		
    end
    
	methods(Static)
		function K = make_kernel_matrix(omega, v, v_b, delta, D, sigma_D)
			% Kernel matrix of the single layer potential L^2(\p D) -> L^2(\p D).
            % Case of two identical boundaries. Although the
            % operator is stil well defined (no jumps), the diagonal
            % terms of the kernel matrix are computed in a special way.
            R_b = D.axis_a;
			green = @(x1,x2,y) tools.G(x1,x2,y,omega,v,v_b,delta,R_b);

            N = size(D,2);
            K = zeros(N) ;

            for i=1:N
                K(1:(i-1),i) = sigma_D(i)*green(D(1,1:(i-1)),D(2,1:(i-1)),[D(1,i),D(2,i)]);
                K((i+1):end,i) = sigma_D(i)*green(D(1,(i+1):end),D(2,(i+1):end),[D(1,i),D(2,i)]);
                
                % Diagonal terms
                K(i,i) = 1/2/pi*sigma_D(i)*(log(sigma_D(i)) - 1);
                %K(i,i) = 1/2 * (K(i, i-1) +  K(i, i+1));
            end
        end
	end
	
end

