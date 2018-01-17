function out = Q_test(k,Num_m)
out = 0;
for a = 1:Num_m
    out = out + 4*besselh(0,1,k*a);
    out = out + 4*besselh(0,1,k*sqrt(2)*a);
    for b = 1:a-1
        absm = sqrt(a^2+b^2);
        out = out+8*besselh(0,1,k*absm);
    end
end