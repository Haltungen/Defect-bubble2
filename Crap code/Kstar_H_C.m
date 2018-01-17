classdef Kstar_H_C < ops.Operators
	% Adjoint of the Neumann-Poincare operator for Helmholtz equation
	
	methods
		function obj = Kstar_H(k, D1, type1, step1)
            
			obj = obj@ops.Operators(D1, type1, step1);
            obj.Kmat = ops.Kstar_H.make_kernel_matrix(k, D1.points, D1.tvec, D1.normal, D1.avec, D1.sigma);

			
		end
		
	end
	
	methods(Static)
		function KsH = make_kernel_matrix(k, D, tvec, normal, avec, sigma)
			% The discretization of the adjoint operator of K using P0 boundary
			% element. Precisely, we construct a matrix Ks such that Ks*phi
			% approximates the integral definition of K^*[phi].
			green = @(x1,x2,y) tools.G(x1,x2,y,omega,v,v_b,delta,R_b);
			M = size(D,2);
			KsH = zeros(M, M);
			tvec_norm_square = tvec(1,:).^2 + tvec(2,:).^2;
			
			for j = 1:M %CHECK!!!!!!!!!!!!!!!!
				xdoty = (D(1,j)-D(1,:))*normal(1,j)+(D(2,j)-D(2,:))*normal(2,j);
				norm_xy = sqrt((D(1,j)-D(1,:)).^2+(D(2,j)-D(2,:)).^2);
                KsH(j,1:j-1) = 1i/4*k*besselh(1,1,k*norm_xy(1:(j-1))).*sigma(1:(j-1)).* xdoty(1:(j-1))./norm_xy(1:(j-1));
				KsH(j,j+1:M) = 1i/4*k*besselh(1,1,k*norm_xy((j+1):M)).*sigma((j+1):M).* xdoty((j+1):M)./norm_xy((j+1):M);
				KsH(j,j) = -1/4/pi*avec(:,j)'*normal(:,j)/tvec_norm_square(j)*sigma(j);

			end
        end
		
		function val = eval()
			error('Method not implemented!');
		end
		
	end
end

