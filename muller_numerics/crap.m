function out = crap(r, rp,n)
if rp > r+0.1
    out = besselh(-n,rp)*besselj(n,r);
else
    out = besselh(n,r)*besselj(-n,rp);
end