function y = E1(x)

% The function E1 as described in [1, Sec. 5.1.53/56]
%
% E1(x) = \int_1^\infty exp(-x*t)/t dt.
%
% [1] A. Abramowitz, I.A. Stegun, Handbook of Mathematical Functions (1970)

%assert(x > 0, 'Evaluation of E1 only for x > 0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
% For x \le 1

A0 = -0.57721566;
A1 =  0.99999193;
A2 = -0.24991055;
A3 =  0.05519968;
A4 = -0.00976004;
A5 =  0.00107857;

%%%%%%%%%%%%%%%%%%%
% For x \ge 1

a1 =  8.5733287401;
a2 = 18.0590169730;
a3 =  8.6347608925;
a4 =  0.2677737343;
b1 =  9.5733223454;
b2 = 25.6329561486;
b3 = 21.0996530827;
b4 =  3.9584969228;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = zeros(size(x));
indLess1 = find(x <= 1);
xL1 = x(indLess1);
y(indLess1) = A0 + A1*xL1 + A2*xL1.^2 + A3*xL1.^3 + A4*xL1.^4 + A5*xL1.^5 - log(xL1);

indGreat1 = find(x > 1);
xG1 = x(indGreat1);

y(indGreat1) = (xG1.^4 + a1*xG1.^3 + a2*xG1.^2 + a3*xG1 + a4)./...
    (xG1.*exp(xG1).*(xG1.^4 + b1*xG1.^3 + b2*xG1.^2 + b3*xG1 + b4));
end
