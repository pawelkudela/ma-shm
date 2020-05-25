function N=shape2Dp(nx,ny,Q,x,y)
 % Calculate shape functions
% nx - number of nodes in 1D shape function in x direction
% ny - number of nodes in 1D shape function in y direction
% Q - Vandermonde 2D matrix
% calculate values of the shape function at point (x,y) defined in local
% coordinates

f=zeros(nx*ny,length(x));
c=0;
for j1=1:ny % eta degree
    for j2=1:nx % ksi degree
        c=c+1;
        f(c,:)=x.^(j2-1).*y.^(j1-1);
    end
end

N=Q'*f;