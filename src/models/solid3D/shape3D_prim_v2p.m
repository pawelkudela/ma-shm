function [Nx,Ny,Nz]=shape3D_prim_v2p(nx,ny,nz,Qx,Qy,Qz,x,y,z)

% nx - number of nodes in 1D shape function in x direction
% ny - number of nodes in 1D shape function in y direction
% nz - number of nodes in 1D shape function in z direction

% calculate values of the shape function at point (x,y,z) defined in local
% coordinates
Nx=zeros(nx*ny*nz,nx*ny*nz);
Ny=zeros(nx*ny*nz,nx*ny*nz);
Nz=zeros(nx*ny*nz,nx*ny*nz);
ux=zeros(nx,nx);
uy=zeros(ny,ny);
uz=zeros(nz,nz);

ux(:,1)=1;
uy(:,1)=1;
uz(:,1)=1;

for j=2:nx
        ux(:,j)=x.^(j-1);
 end;  
 for j=2:ny
        uy(:,j)=y.^(j-1);
 end
 for j=2:nz
        uz(:,j)=z.^(j-1);
 end;  
 Sx=ux*Qx';
 Sy=uy*Qy';
 Sz=uz*Qz';
 
uxprim=zeros(nx,nx);
uyprim=zeros(ny,ny);
uzprim=zeros(nz,nz);

for j=2:nx
        uxprim(:,j)=(j-1)*x.^(j-2);
 end;  
 for j=2:ny
        uyprim(:,j)=(j-1)*y.^(j-2);
 end
 for j=2:nz
        uzprim(:,j)=(j-1)*z.^(j-2);
 end;  
 Sxprim=uxprim*Qx';
 Syprim=uyprim*Qy';
 Szprim=uzprim*Qz';
 c=0;
 for k3=1:nz
 for k2=1:ny
 for k1=1:nx
 cc=0;c=c+1;
 for j3=1:nz
     for j2=1:ny % shape function number
         for j1=1:nx
             cc=cc+1;                 
             Nx(c,cc)=Sxprim(k1,j1)*Sy(k2,j2)*Sz(k3,j3);
             if(abs(Nx(c,cc))<1e-12) Nx(c,cc)=0; end
             Ny(c,cc)=Sx(k1,j1)*Syprim(k2,j2)*Sz(k3,j3);
             if(abs(Ny(c,cc))<1e-12) Ny(c,cc)=0; end
             Nz(c,cc)=Sx(k1,j1)*Sy(k2,j2)*Szprim(k3,j3);
             if(abs(Nz(c,cc))<1e-12) Nz(c,cc)=0; end
          end
     end
 end
 end
 end
 end
         
   

