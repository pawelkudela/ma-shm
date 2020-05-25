function N=shape3D_v2(nx,ny,nz,Qx,Qy,Qz,x,y,z)

% nx - number of nodes in 1D shape function in x direction
% ny - number of nodes in 1D shape function in y direction
% nz - number of nodes in 1D shape function in z direction

% calculate values of the shape function at point (x,y,z) defined in local
% coordinates

N=zeros(1,nx*ny*nz);
ux=zeros(1,nx);
uy=zeros(1,ny);
uz=zeros(1,nz);
ux(1)=1;
uy(1)=1;
uz(1)=1;

for j=2:nx
        ux(j)=x^(j-1);
 end;  
 for j=2:ny
        uy(j)=y^(j-1);
 end
 for j=2:nz
        uz(j)=z^(j-1);
 end;  
 Sx=ux*Qx';
 Sy=uy*Qy';
 Sz=uz*Qz';
 cc=0;
 for j3=1:nz
     for j2=1:ny % shape function number
         for j1=1:nx
             cc=cc+1;                 
             N(cc)=Sx(j1)*Sy(j2)*Sz(j3);
          end
     end
end
         
   

