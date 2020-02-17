function [kufi,kfifi]=coupl3D_global(epzt,gpzt,x,y,z,Qx,Qy,Qz,ksi,eta,dzeta,wx,wy,wz)

%calculate electrical matrices for solid spectral element



nx=length(ksi);
ny=length(eta);
nz=length(dzeta);

% 3D case 
B=zeros(6,nx*ny*nz*3);
Bfi=zeros(3,nx*ny*nz);
kufi=zeros(nx*ny*nz*3,nx*ny*nz);
kfifi=zeros(nx*ny*nz,nx*ny*nz);
for kz=1:nz % dzeta
    cm=0;
    for ky=1:ny % eta
     for kx=1:nx % ksi     
         J=zeros(3,3);
         kk=0;
         [Nprimx,Nprimy,Nprimz]=shape3D_prim_v2(nx,ny,nz,Qx,Qy,Qz,ksi(kx),eta(ky),dzeta(kz));
         %         S=shape3D_v2(nx,ny,nz,Qx,Qy,Qz,ksi,eta,dzeta);
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
%          if(abs(J(1,1)<eps)) J(1,1)=0; end;
%          if(abs(J(1,2)<eps)) J(1,2)=0; end;
%          if(abs(J(1,3)<eps)) J(1,3)=0; end;
%          if(abs(J(2,1)<eps)) J(2,1)=0; end;
%          if(abs(J(2,2)<eps)) J(2,2)=0; end;
%          if(abs(J(2,3)<eps)) J(2,3)=0; end;
%          if(abs(J(3,1)<eps)) J(3,1)=0; end;
%          if(abs(J(3,2)<eps)) J(3,2)=0; end;
%          if(abs(J(3,3)<eps)) J(3,3)=0; end;
         Jj=[J(1,1),J(1,2),J(1,3),J(2,1),J(2,2),J(2,3),J(3,1),J(3,2),J(3,3)];
         
         % normal vector
         V3(1)=Jj(2)*Jj(6)-Jj(5)*Jj(3); 
         V3(2)=Jj(4)*Jj(3)-Jj(1)*Jj(6);
         V3(3)=Jj(1)*Jj(5)-Jj(4)*Jj(2);  
         y0=V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3);
         ys=sqrt(y0);
         v31=V3(1)/ys;
         v32=V3(2)/ys;
         v33=V3(3)/ys;
         %/*v1 vector*/
         V1(1)=Jj(1);
         V1(2)=Jj(2);
         V1(3)=Jj(3);
         y0=V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3); 
         ys=sqrt(y0);
         v11=V1(1)/ys;
         v12=V1(2)/ys;
         v13=V1(3)/ys;
         %/*cross product v2=v3xv1*/
         v21=v32*v13-v33*v12;
         v22=-v31*v13+v33*v11;
         v23=v31*v12-v32*v11;
         y0=v21*v21+v22*v22+v23*v23;
         ys=sqrt(y0);
         v21=v21/ys;
         v22=v22/ys;
         v23=v23/ys;
         %%%
         % /* determinant and inverse Jacobian  */
         inJ=inv(J);
         detJ=det(J);
         wwwdetJ=wx(kx)*wy(ky)*wz(kz)*detJ;
         cc=0;
         for j3=1:nz
          for j2=1:ny % shape function number
              for j1=1:nx
                  cc=cc+1;                 
                Nx=Nprimx(cc)*inJ(1,1)+Nprimy(cc)*inJ(1,2)+Nprimz(cc)*inJ(1,3);
                Ny=Nprimx(cc)*inJ(2,1)+Nprimy(cc)*inJ(2,2)+Nprimz(cc)*inJ(2,3);
                Nz=Nprimx(cc)*inJ(3,1)+Nprimy(cc)*inJ(3,2)+Nprimz(cc)*inJ(3,3);
                %B - strain - nodal displacement matrix;
                B11=Nx;
                B12=0;
                B13=0;
                B21=0;
                B22=Ny;
                B23=0;
                B31=0;
                B32=0;
                B33=Nz;
                B41=Ny;
                B42=Nx;
                B43=0;
                B51=0;
                B52=Nz;
                B53=Ny;
                B61=Nz;
                B62=0;
                B63=Nx;

                % epsx=du/dx
                B(1,3*(cc)-2)=B11;
                B(1,3*(cc)-1)=B12;
                B(1,3*(cc))=B13;
                % epsy=dv/dy
                B(2,3*(cc)-2)=B21;
                B(2,3*(cc)-1)=B22;
                B(2,3*(cc))=B23;
                % epsz=dw/dz
                B(3,3*(cc)-2)=B31; 
                B(3,3*(cc)-1)=B32;
                B(3,3*(cc))=B33;
                % gammaxy=du/dy+dv/dx
                B(4,3*(cc)-2)=B41; 
                B(4,3*(cc)-1)=B42;
                B(4,3*(cc))=B43;
                % gammayz=dw/dy+dv/dz
                B(5,3*(cc)-2)=B51; 
                B(5,3*(cc)-1)=B52;
                B(5,3*(cc))=B53;
                % gammaxz=dw/dx+du/dz
                B(6,3*(cc)-2)=B61; 
                B(6,3*(cc)-1)=B62;
                B(6,3*(cc))=B63;
                %%
                Bfi(1,cc)=Nx*v11 + Ny*v12 + Nz*v13;
                Bfi(2,cc)=Nx*v21 + Ny*v22 + Nz*v23;
                Bfi(3,cc)=Nx*v31 + Ny*v32 + Nz*v33;
              end
          end
         end     
         cm=cm+1;
         kufi=kufi+(B'*epzt'*Bfi)*wwwdetJ;
         kfifi=kfifi+Bfi'*gpzt*Bfi*wwwdetJ;
     end
    end
end

