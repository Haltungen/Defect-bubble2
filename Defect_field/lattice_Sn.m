function Sn=lattice_Sn(n,k,alp,L,Num_m)


    
    
if mod(n,2) == 0
    %number is even
      
    if n==0
      Sn=S1D_0(k,alp,L,Num_m)+DeltaSn(0,k,alp,L,L,Num_m);  
    else
      l=n/2;
      Sn=S1D_2n(l,k,alp,L,Num_m)+DeltaSn(n,k,alp,L,L,Num_m);
      %Sn=Sn*(-1)^n;
    end
else
  %number is odd
  l=(n+1)/2;
  Sn=S1D_2n_1(l,k,alp,L,Num_m)+DeltaSn(n,k,alp,L,L,Num_m);
  %Sn=Sn*(-1)^n;
end




