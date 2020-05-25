function N=shape1D_v2(nz,Qz,z)

% nz - number of nodes in 1D shape function in z direction

% calculate values of the shape function at point (z) defined in local
% coordinates

N=zeros(1,nz);
uz=zeros(1,nz);
uz(1)=1;

 for j=2:nz
        uz(j)=z^(j-1);
 end;  
 Sz=uz*Qz';
 cc=0;
 for j3=1:nz
     cc=cc+1;                 
     N(cc)=Sz(j3);
end
         
   

