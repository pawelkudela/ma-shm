close all;clear all;

for k_test=[870];

%Nx=400; %Nx=152;% grid size
%Ny=400; %Ny=300;
% 29
Nx=500;%
Ny=500;
tic

Eps=1e-8;
mode='gpu';% options: mode='cpu';mode='gpu';
% currently only for bottom surface
surface='bottom'; % options surface='top';surface='bottom';surface='middle';

os='win';
switch os
    case 'win'
        input_file = ['inputs\input',num2str(k_test)];
        run(input_file); 
        data_filename=['outputs\output',num2str(k_test),'\solid3D_in_plane_',num2str(k_test),'_',num2str(Nx),'x',num2str(Ny),surface,'.mat'];
    case 'linux'
        input_file =['inputs/input',num2str(k_test)]; 
        run(input_file);
        data_filename=['outputs/output',num2str(k_test),'/solid3D_in_plane_',num2str(k_test),'_',num2str(Nx),'x',num2str(Ny),surface,'.mat'];

    otherwise
        disp('Unknown system');
end
%nft=nft/2; % for testing (half of the time)
%frm_int=2*frm_int;
switch mode
    case 'gpu'
        Data=zeros(Ny,Nx,nft/frm_int,'double','gpuArray');
    case 'cpu'
        Data=zeros(Ny,Nx,nft/frm_int);
end
load(meshfile);
switch mode
    case 'gpu'
        %nodes=gpuArray(uint32(nodes));
        nodes=gpuArray(nodes);
        coords=gpuArray(coords);
    case 'cpu'
        
end
[fen,no]=size(nodes);
lay=1;
%fen=fen-2*length(pztnum);
fen=fen-1*length(pztnum)-length(stiffenerEl);% no glue
L=max(coords(:,1))-min(coords(:,1))-2*Eps;
xi=min(coords(:,1))+Eps:L/(Nx-1):max(coords(:,1))-Eps;
B=max(coords(:,2))-min(coords(:,2))-2*Eps;
yi=min(coords(:,2))+Eps:B/(Ny-1):max(coords(:,2))-Eps;
[XI,YI]=meshgrid(xi,yi);
xp=reshape(XI,1,Nx*Ny);
yp=reshape(YI,1,Nx*Ny);
switch mode
    case 'gpu'
        ZA=zeros(nx*ny,Nx*Ny,'double','gpuArray');
        %ZI=zeros(Nx,Ny,'double','gpuArray');
    case 'cpu'
        ZA=zeros(nx*ny,Nx*Ny);
        %ZI=zeros(Nx,Ny);
end


switch surface
   case 'bottom'
   ne=1:fen/lay;
   case 'top'
   ne=(fen/lay)*(lay-1)+1:fen;
   case 'middle'
   ne=(fen/lay)*(lay-3)+1:fen;  
end
disp('Preparing shape functions');
%%
% extract corner nodes of all elements
x=coords(nodes(ne,[1,6,36,31]),1);
y=coords(nodes(ne,[1,6,36,31]),2);

x=reshape(x,fen,4)';
y=reshape(y,fen,4)';
 
% find centres (xc,yc) of all elements
% bx=(x(4,:)+x(3,:))/2;
% by=(y(4,:)+y(3,:))/2;
% ax=(x(3,:)+x(2,:))/2;
% ay=(y(3,:)+y(2,:))/2;
% 
% cx=(x(1,:)+x(4,:))/2;
% cy=(y(1,:)+y(4,:))/2;
% dx=(x(2,:)+x(1,:))/2;
% dy=(y(2,:)+y(1,:))/2;
% 
% xc = (ax+cx)/2;
% yc = (by+dy)/2;

% cross product method
% check if there is the same sign so that the point is inside convex polygon (element)
switch mode
    case 'gpu'
        elem_index=zeros(Nx*Ny,1,'uint32','gpuArray');
        %elem_index=zeros(Nx*Ny,1,'double','gpuArray');
    case 'cpu'
        elem_index=zeros(Nx*Ny,1);
end

x21=x(2,:)-x(1,:);
x32=x(3,:)-x(2,:);
x43=x(4,:)-x(3,:);
x14=x(1,:)-x(4,:);
y21=y(2,:)-y(1,:);
y32=y(3,:)-y(2,:);
y43=y(4,:)-y(3,:);
y14=y(1,:)-y(4,:);
x1y21=x(1,:).*y21;
y1x21=y(1,:).*x21;
x2y32=x(2,:).*y32;
y2x32=y(2,:).*x32;
x3y43=x(3,:).*y43;
y3x43=y(3,:).*x43;
x4y14=x(4,:).*y14;
y4x14=y(4,:).*x14;
disp('Find elements');
switch mode
    case 'gpu'
        poolobj = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(poolobj)
            poolobj = parpool; % create parpool
        end
        parfor j =1:Nx*Ny
        %for j =1:Nx*Ny
            % cross product for the consequtive edge
            d1 = (xp(j).*y21-x1y21)  - (yp(j).*x21-y1x21);
            d2 = (xp(j).*y32-x2y32)  - (yp(j).*x32-y2x32);
            d3 = (xp(j).*y43-x3y43)  - (yp(j).*x43-y3x43);
            d4 = (xp(j).*y14-x4y14)  - (yp(j).*x14-y4x14);
    
            negd1=uint32(d1<0);negd2=uint32(d2<0);negd3=uint32(d3<0);negd4=uint32(d4<0);
            %negd1=(d1<0);negd2=(d2<0);negd3=(d3<0);negd4=(d4<0);
            neg=negd1.*negd2.*negd3.*negd4;
            neg_index=find(neg);
            elem_index(j,1)=neg_index(1);
        end
        delete(poolobj);
    case 'cpu'
        for j =1:Nx*Ny
            % cross product for the consequtive edge
            d1 = (xp(j).*y21-x1y21)  - (yp(j).*x21-y1x21);
            d2 = (xp(j).*y32-x2y32)  - (yp(j).*x32-y2x32);
            d3 = (xp(j).*y43-x3y43)  - (yp(j).*x43-y3x43);
            d4 = (xp(j).*y14-x4y14)  - (yp(j).*x14-y4x14);
    
            negd1=uint32(d1<0);negd2=uint32(d2<0);negd3=uint32(d3<0);negd4=uint32(d4<0);
            neg=negd1.*negd2.*negd3.*negd4;
            neg_index=find(neg);
            elem_index(j,1)=neg_index(1);
        end
end
       
disp('Calculate local values of ksi and eta for each point');
% solve equation system: XI = N1*x1+N2*x2+N3*x3+N4*x4
%                                 YI = N1*y1+N2*y2+N3*y3+N4*y4
% in which (x,y) are corner nodes in global coordinate system
% (XI,YI) are coordinates of selected points (on regular grid)
% N are shape functions of quad element
% N1 = 1/4*(1-ksi)*(1-eta)
% N2 = 1/4*(1+ksi)*(1-eta)
% N3 = 1/4*(1+ksi)*(1+eta)
% N4 = 1/4*(1-ksi)*(1+eta)
% ksi, eta are unknowns

sx=x(1,elem_index)+x(2,elem_index)+x(3,elem_index)+x(4,elem_index);
sy=y(1,elem_index)+y(2,elem_index)+y(3,elem_index)+y(4,elem_index);
rx=-x(1,elem_index)-x(2,elem_index)+x(3,elem_index)+x(4,elem_index);
ry=-y(1,elem_index)-y(2,elem_index)+y(3,elem_index)+y(4,elem_index);
px=-x(1,elem_index)+x(2,elem_index)+x(3,elem_index)-x(4,elem_index);
py=-y(1,elem_index)+y(2,elem_index)+y(3,elem_index)-y(4,elem_index);
qx=x(1,elem_index)-x(2,elem_index)+x(3,elem_index)-x(4,elem_index);
qy=y(1,elem_index)-y(2,elem_index)+y(3,elem_index)-y(4,elem_index);

A1=qy.*rx-ry.*qx;
B1=4*yp.*qx-sy.*qx-ry.*px+py.*rx+sx.*qy-4*xp.*qy;
C1=4*yp.*px-sy.*px-4*py.*xp+py.*sx;
switch mode
    case 'gpu'
        ksi_p1 = zeros(1,Nx*Ny,'double','gpuArray');
        eta_p1= zeros(1,Nx*Ny,'double','gpuArray');
        ksi_p2 = zeros(1,Nx*Ny,'double','gpuArray');
        eta_p2= zeros(1,Nx*Ny,'double','gpuArray');
        ksi_p = zeros(1,Nx*Ny,'double','gpuArray');
        eta_p= zeros(1,Nx*Ny,'double','gpuArray');
        %elem_index=zeros(Nx*Ny,1,'double','gpuArray');
    case 'cpu'
        ksi_p1 = zeros(1,Nx*Ny);
        eta_p1= zeros(1,Nx*Ny);
        ksi_p2 = zeros(1,Nx*Ny);
        eta_p2= zeros(1,Nx*Ny);
        ksi_p = zeros(1,Nx*Ny);
        eta_p= zeros(1,Nx*Ny);
end

II=find(A1==0); % linear equation
eta_p1(II) = -C1(II)./B1(II);
ksi_p1(II)=( 4*xp(II)-sx(II)-eta_p1(II).*rx(II) )./(px(II)+eta_p1(II).*qx(II));
% resolve infinite ksi
Iinf=find(isinf(ksi_p1));
ksi_p1(Iinf)=( 4*yp(Iinf) -sy(Iinf)-eta_p1(Iinf).*ry(Iinf) )./( py(Iinf) + eta_p1(Iinf).*qy(Iinf));
% resolve nan ksi
Inan=find(isnan(ksi_p1));
ksi_p1(Inan) = 0;

JJ=setdiff([1:Nx*Ny],II); % quadratic equation
Delta1=B1(JJ).^2-4*A1(JJ).*C1(JJ);
eta_p1(JJ)= ( -B1(JJ) -sqrt(Delta1) )./(2*A1(JJ));
eta_p2(JJ)= ( -B1(JJ) +sqrt(Delta1) )./(2*A1(JJ));
ksi_p1(JJ)=( 4*xp(JJ)-sx(JJ)-eta_p1(JJ).*rx(JJ) )./(px(JJ)+eta_p1(JJ).*qx(JJ));
ksi_p2(JJ)=( 4*xp(JJ)-sx(JJ)-eta_p2(JJ).*rx(JJ) )./(px(JJ)+eta_p2(JJ).*qx(JJ));
% select only appropriate solutions of quadratic equation
for k=1:Nx*Ny
    if( abs(eta_p1(k))<=1+Eps )
        eta_p(k) = eta_p1(k);
        ksi_p(k) = ksi_p1(k);
    else
        eta_p(k) = eta_p2(k);
        ksi_p(k) = ksi_p2(k);
    end
end

% evaluate shape functions
disp('Evaluate shape functions');
switch mode
    case 'gpu'
        [N]=shape_func_p_gpu(ksi_p,eta_p);
    case 'cpu'
        [N]=shape_func_p(ksi_p,eta_p);
end
clear ksi_p eta_p ksi_p1 eta_p1 ksi_p2 eta_p2
% create sparse matrix
% Ns1 = sparse(Nx*Ny,nx*ny*Nx*Ny);
% for k=1:Nx*Ny
%     Ns1(k,(k-1)*nx*ny+1:k*nx*ny) = N(:,k)';
% end
%
[indxi,indxj]=find(N);
Ns=sparse(indxj,[1:nx*ny*Nx*Ny],reshape(N,nx*ny*Nx*Ny,1));
clear N;
%%
sample=0;
for n=frm_int:frm_int:nft
    sample=sample+1;
    [sample,nft/frm_int]
switch os
    case 'win'
        Ux_filename=['outputs\output',num2str(k_test),'\Ux_frame',num2str(n),'.mat'];
        Uy_filename=['outputs\output',num2str(k_test),'\Uy_frame',num2str(n),'.mat'];
    case 'linux'
        Ux_filename=['outputs/output',num2str(k_test),'/Ux_frame',num2str(n),'.mat'];
        Uy_filename=['outputs/output',num2str(k_test),'/Uy_frame',num2str(n),'.mat'];
        
    otherwise
        disp('Unknown system');
end


load(Ux_filename,'Ux'); % fastest
load(Uy_filename,'Uy'); % fastest

switch mode
    case 'gpu'
        Ux=gpuArray(Ux);
        Uy=gpuArray(Uy);
    case 'cpu'
        
end

%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch surface
   case 'bottom'
      UX=Ux(nodes(elem_index,1:nx*ny))'; % bottom surface
      UY=Uy(nodes(elem_index,1:nx*ny))'; % bottom surface
   case 'top'
      UX=Ux(nodes(elem_index,nx*ny*nz-nx*ny+1:nx*ny*nz))'; % top surface
      UY=Uy(nodes(elem_index,nx*ny*nz-nx*ny+1:nx*ny*nz))'; % top surface
   case 'middle'
      UX=Ux(nodes(elem_index,nx*ny*nz-nx*ny+1:nx*ny*nz))'; % middle surface
      UY=Uy(nodes(elem_index,nx*ny*nz-nx*ny+1:nx*ny*nz))'; % middle surface
end
% c=0;
% for j =1:Ny
%     for i=1:Nx
%         c=c+1;
%         ZI(i,j) = N(:,c)'*ZA(:,c);
%     end
% end
ZA=sqrt(UX.^2+UY.^2);
ZA=reshape(ZA,nx*ny*Nx*Ny,1);
ZI=Ns*ZA;% interpolated values
ZI=reshape(ZI,Nx,Ny);
Data(:,:,sample)=ZI;

end

switch mode
    case 'gpu'
        Data=gather(Data);
        
    case 'cpu'
        
end

toc
save(data_filename,'Data');
%save(data_filename,'Data','-v7.3'); % for data files >2GB (compression is used)
pause(0.1);
end