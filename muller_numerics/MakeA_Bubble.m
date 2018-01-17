function A = MakeA_Bubble(omega,v0,vb,alp,delta,L,R,NN,N)
k0 = omega*v0;
kb = omega*vb;

M=2*NN+1;
A11=zeros(M,M);
A21=zeros(M,M);
A12c=zeros(M,M);
A12d=zeros(M,M);
A22c=zeros(M,M);
A22d=zeros(M,M);

%mm=1;


%%% store lattice sum data %%%
data_Sn_posi = zeros(2*NN+1,1);
data_Sn_nega = zeros(2*NN,1);
for j=1:2*NN+1
data_Sn_posi(j)=tools.lattice_Sn(j-1,k0,alp,L,N);
end
for j=1:2*NN   
data_Sn_nega(j)=tools.lattice_Sn(j,k0,[-alp(1), alp(2)],L,N);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=-NN:NN
    J0 = besselj(n, k0*R);
    H0 = besselh(n, 1, k0*R);
    
    Jb = besselj(n, kb*R);
    Hb = besselh(n, 1, kb*R);
    
    dH0 = 1/2*(besselh(n-1,1, k0*R) - besselh(n+1,1, k0*R));
    dJb = 1/2*(besselj(n-1, kb*R) - besselj(n+1, kb*R));

    c=-0.5*1i*pi*R;
    
    %inside
    A11(NN+1+n,NN+1+n)=c*Jb*Hb;
    A21(NN+1+n,NN+1+n)=c*kb*dJb*Hb;
    
    %outside
    A12c(NN+1+n,NN+1+n)=c*J0*H0;
    A22c(NN+1+n,NN+1+n)=c*k0*J0*dH0;
    
    for m=-NN:NN
        if (n-m) >= 0
            Snm0 = data_Sn_posi(n-m+1);
        else
            Snm0 = data_Sn_nega(-(n-m));
        end
        dJ0m = 1/2*(besselj(m-1,k0*R) - besselj(m+1,k0*R));
        J0m = besselj(m,k0*R);
        
        A12d(NN+1+m,NN+1+n)=c*J0*(-1)^(n-m)*Snm0*J0m;
        A22d(NN+1+m,NN+1+n)=c*k0*J0*(-1)^(n-m)*Snm0*dJ0m;
    end
end

A=[A11,-A12c-A12d; A21, delta*(A22c+A22d)];
end



