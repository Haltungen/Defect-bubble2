function sum=DeltaSn(n,k,alp,L,L2,Num_m)



sum=0;




for m=-Num_m:Num_m
    
    
    bm=alp(1)+m*2*pi/L;
   % bm2=alp(1)+(-m)*2*pi/L;
    
    tm=asin(bm/k);
    %tm2=asin(bm2/k);

    
    gm=sqrt(k^2-bm^2);
    %gm2=sqrt(bm2^2-k^2);
    
    
    sum=sum+1/gm*(  exp(1i*n*tm)/(exp(-1i*alp(2)*L2)*exp(-1i*gm*L2) - 1)  + (-1)^n*exp(-1i*n*tm)/(exp(1i*alp(2)*L2)*exp(-1i*gm*L2) - 1)  );
    
    
    
    

end

sum=sum*2/L;


