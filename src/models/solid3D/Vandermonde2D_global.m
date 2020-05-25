function [Q]=Vandermonde2D_global(X,Y,nx,ny)
% Calculate coefficients for Vandermonde matrices for shape function representation
% n - number of nodes


V=zeros(nx*ny,nx*ny);
r=0;
for i1=1:ny % eta point
    for i2=1:nx%ksi point
        c=0;r=r+1;
        for j1=1:ny % eta degree
            for j2=1:nx % ksi degree
                c=c+1;
                V(r,c)=X(r)^(j2-1)*Y(r)^(j1-1);
            end
        end
    end
end
Q=inv(V);
