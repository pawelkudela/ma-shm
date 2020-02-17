function [Vx,Vzx]=Vandermonde_old(ksi,ksiz,n,nz)

V=zeros(n,n);
for i=1:n
    for j=1:n
        V(i,j)=ksi(j)^(i-1);
    end;
end
Q=inv(V);
Vx=zeros(n,n);
for k=1:n
    for i=1:n
        for j=2:n
            Vx(k,i)=Vx(k,i)+Q(i,j)*(j-1)*ksi(k)^(j-2);
        end;
    end
end
V=zeros(n,n);
for k=1:n
    V(k,k)=1;
end;
for k=2:n-1
    Vx(k,k)=0;
end

% Gauss-Lobatto_Legendre points and weights
% Vandermonde matrix
Vz=zeros(nz,nz);
for i=1:nz
    for j=1:nz
        Vz(i,j)=ksiz(j)^(i-1);
    end;
end
Qz=inv(Vz);
Vzx=zeros(nz,nz);
for k=1:nz
    for i=1:nz
        for j=2:nz
            Vzx(k,i)=Vzx(k,i)+Qz(i,j)*(j-1)*ksiz(k)^(j-2);
        end;
    end
end
Vz=zeros(nz,nz);
for k=1:nz
    Vz(k,k)=1;
end;
for k=2:nz-1
    Vzx(k,k)=0;
end