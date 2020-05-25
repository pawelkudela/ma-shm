function [Nx,Ny]=shape_derivatives2Dp(nx,ny,Q,x,y)
 % Calculate shape function derivatives
% nx - number of nodes in 1D shape function in x direction
% ny - number of nodes in 1D shape function in y direction
% Q - Vandermonde 2D matrix
% Nx - derivative in respect to x
% Ny - derivative in respect to y
% calculate values of the shape function derivatives at point (x,y) defined in local
% coordinates


fx=zeros(nx*ny,length(x));
c=0;
for j1=1:ny % eta degree
    for j2=1:nx % ksi degree
        c=c+1;
        fx(c,:)=((j2-1)*x.^(j2-2)).*y.^(j1-1);
    end
end

Nx=Q'*fx;

fy=zeros(nx*ny,length(x));
c=0;
for j1=1:ny % eta degree
    for j2=1:nx % ksi degree
        c=c+1;
        fy(c,:)=x.^(j2-1).*((j1-1)*y.^(j1-2));
    end
end

Ny=Q'*fy;