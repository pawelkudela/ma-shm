function N=shape2D(nx,ny,Qx,Qy,x,y)
 % Calculate shape functions at GLL nodes
% nx - number of nodes in 1D shape function in x direction
% ny - number of nodes in 1D shape function in y direction

% calculate values of the shape function at point (x,y) defined in local
% coordinates

N=zeros(nx*ny,nx*ny);

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

 c=0;
 
 for k2=1:ny
 for k1=1:nx
 cc=0;c=c+1;
 
     for j2=1:ny % shape function number
         for j1=1:nx
             cc=cc+1;                 
             N(c,cc)=Sx(k1,j1)*Sy(k2,j2);
          end
     end

 end
 end

         
   

