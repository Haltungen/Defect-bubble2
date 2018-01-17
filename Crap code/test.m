h = 0.01;
func1 = @(x) (besselh(0,x*(1+h))-2*besselh(0,x)+besselh(0,x*(1-h)))./(x*h).^2;
func2 = @(x) -besselh(0,x)+1./x.*besselh(1,x);
x = linspace(0.1,1);
plot(x,func1(x),'*',x,func2(x))