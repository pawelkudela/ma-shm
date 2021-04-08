function [Nx,Ny,Nz]=shape3D_prim_v2(nx,ny,nz,Qx,Qy,Qz,x,y,z)

% nx - number of nodes in 1D shape function in x direction
% ny - number of nodes in 1D shape function in y direction
% nz - number of nodes in 1D shape function in z direction

% calculate values of the shape function at point (x,y,z) defined in local
% coordinates
Nx=zeros(1,nx*ny*nz);
Ny=zeros(1,nx*ny*nz);
Nz=zeros(1,nx*ny*nz);
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
 
uxprim=zeros(1,nx);
uyprim=zeros(1,ny);
uzprim=zeros(1,nz);

for j=2:nx
        uxprim(j)=(j-1)*x^(j-2);
 end;  
 for j=2:ny
        uyprim(j)=(j-1)*y^(j-2);
 end
 for j=2:nz
        uzprim(j)=(j-1)*z^(j-2);
 end;  
 Sxprim=Qx*uxprim';
 Syprim=uyprim*Qy';
 Szprim=uzprim*Qz';
 cc=0;
 for j3=1:nz
     for j2=1:ny % shape function number
         for j1=1:nx
             cc=cc+1;                 
             Nx(cc)=Sxprim(j1)*Sy(j2)*Sz(j3);
             if(abs(Nx(cc))<1e-12) Nx(cc)=0; end
             Ny(cc)=Sx(j1)*Syprim(j2)*Sz(j3);
             if(abs(Ny(cc))<1e-12) Ny(cc)=0; end
             Nz(cc)=Sx(j1)*Sy(j2)*Szprim(j3);
             if(abs(Nz(cc))<1e-12) Nz(cc)=0; end
          end
     end
end
         
   

