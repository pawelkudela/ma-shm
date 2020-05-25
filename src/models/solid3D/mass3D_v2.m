function mg=mass3D_v2(rho,x,y,z,Qx,Qy,Qz,ksi,eta,dzeta,wx,wy,wz)

%calculate mass matrix for solid spectral element
%using Nodal Reduced Quadrature (integration at GLL points)
%mass density is averaged hence rho is scalar


nx=length(ksi);
ny=length(eta);
nz=length(dzeta);

% 3D case 
mg=zeros(nx*ny*nz*3,1);
% mg=sparse(nx*ny*nz*3,nx*ny*nz*3);
% N=sparse(3,nx*ny*nz*3);
c=0;
for kz=1:nz % dzeta
    for ky=1:ny % eta
     for kx=1:nx % ksi     
         J=zeros(3,3);
         kk=0;c=c+1;
         [Nprimx,Nprimy,Nprimz]=shape3D_prim_v2(nx,ny,nz,Qx,Qy,Qz,ksi(kx),eta(ky),dzeta(kz));
         S=shape3D_v2(nx,ny,nz,Qx,Qy,Qz,ksi(kx),eta(ky),dzeta(kz));
         for j3=1:nz
          for j2=1:ny % shape function number
              for j1=1:nx
                  kk=kk+1;
                  J(1,1)=J(1,1)+Nprimx(kk)*x(kk);
                  J(1,2)=J(1,2)+Nprimx(kk)*y(kk);
                  J(1,3)=J(1,3)+Nprimx(kk)*z(kk);
                  
                  J(2,1)=J(2,1)+Nprimy(kk)*x(kk);
                  J(2,2)=J(2,2)+Nprimy(kk)*y(kk);
                  J(2,3)=J(2,3)+Nprimy(kk)*z(kk);
                  
                  J(3,1)=J(3,1)+Nprimz(kk)*x(kk);
                  J(3,2)=J(3,2)+Nprimz(kk)*y(kk);
                  J(3,3)=J(3,3)+Nprimz(kk)*z(kk);
              end
          end
         end
         % /* determinant of Jacobi matrix  */
     
         detJ=det(J);
        
         wwwdetJ=wx(kx)*wy(ky)*wz(kz)*detJ;
%          cc=0;
%          for j3=1:nz
%           for j2=1:ny % shape function number
%               for j1=1:nx
%                   cc=cc+1;                     
%                 % 
%                 N(1,3*(cc)-2)=S(cc);
%                 N(2,3*(cc)-1)=S(cc);
%                 N(3,3*(cc))=S(cc); 
%               end
%           end
%          end     
%          mg=mg+N'*N*wwwdetJ*rho;
           
           mg(3*c-2)=wwwdetJ*rho;
           mg(3*c-1)=wwwdetJ*rho;
           mg(3*c)=wwwdetJ*rho;
     end
    end
end
% mg=diag(mg);
% mg=full(mg);
