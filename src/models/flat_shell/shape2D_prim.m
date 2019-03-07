function [Nx,Ny]=shape2D_prim(nx,ny,Qx,Qy,x,y)
    % Calculate shape function derivatives at GLL nodes
    % nx - number of nodes in 1D shape function in x direction
    % ny - number of nodes in 1D shape function in y direction


% calculate values of the shape function at point (x, y) defined in
% local coordinates (elemental coodinates)

Nx=zeros(nx*ny,nx*ny);
Ny=zeros(nx*ny,nx*ny);

ux=zeros(nx,nx);
uy=zeros(ny,ny);


ux(:,1)=1;
uy(:,1)=1;

for j=2:nx
        ux(:,j)=x.^(j-1);
end  
 for j=2:ny
        uy(:,j)=y.^(j-1);
 end

 Sx=ux*Qx';
 Sy=uy*Qy';
 
uxprim=zeros(nx,nx);
uyprim=zeros(ny,ny);


for j=2:nx
        uxprim(:,j)=(j-1)*x.^(j-2);
end  
 for j=2:ny
        uyprim(:,j)=(j-1)*y.^(j-2);
 end

 Sxprim=uxprim*Qx';
 Syprim=uyprim*Qy';
 c=0;

 for k2=1:ny
 for k1=1:nx
 cc=0;c=c+1;
     for j2=1:ny % shape function number
         for j1=1:nx
             cc=cc+1;                 
             Nx(c,cc)=Sxprim(k1,j1)*Sy(k2,j2);
             if(abs(Nx(c,cc))<1e-12) Nx(c,cc)=0; end
             Ny(c,cc)=Sx(k1,j1)*Syprim(k2,j2);
             if(abs(Ny(c,cc))<1e-12) Ny(c,cc)=0; end
            
          end
     end
 end
 end

         
   

