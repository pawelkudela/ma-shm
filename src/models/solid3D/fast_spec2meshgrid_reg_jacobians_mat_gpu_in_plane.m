close all;clear all;
% only for bottom surface
for k_test=[5:6]

%Nx=400; %Nx=152;% grid size
%Ny=400; %Ny=300;
% 29
Nx=500;%
Ny=500;
% Nx=300;%
% Ny=300;

tic

Eps=1e-8;
mode='cpu';% options: mode='cpu';mode='gpu';
% currently only for bottom surface
surface='bottom'; % options surface='top';surface='bottom';surface='middle';
input_file = fullfile('inputs',['input',num2str(k_test)]);
data_filename=fullfile('outputs',['\output',num2str(k_test)],['\solid3D_in_plane_',num2str(k_test),'_',num2str(Nx),'x',num2str(Ny),surface,'.mat']);
run(input_file); 
load([meshfile(1:end-4),'_jacobians']);
%nft=nft/2; % for testing (half of the time)
%frm_int=2*frm_int;

switch mode
    case 'gpu'
        Data=zeros(Ny,Nx,nft/frm_int,'double','gpuArray');
    case 'cpu'
        Data=zeros(Ny,Nx,nft/frm_int);
end
load(meshfile);
coords=coords3D;
nodes=nodes3D;
NofElNodesx=nx;
NofElNodesy=ny;

switch mode
    case 'gpu'
        %nodes=gpuArray(uint32(nodes));
        nodes=gpuArray(nodes);
        coords=gpuArray(coords);
    case 'cpu'
        
end
[fen,NofElNodes]=size(nodes);
%lay=1;
%fen=fen-2*length(pztnum);
fen=size(msh.QUADS,1);
%fen=fen-1*length(pztnum)-length(stiffenerEl);% no glue
%fen=fen-1*length(pztnum);% no glue
L=max(coords(:,1))-min(coords(:,1))-2*Eps;
xi=min(coords(:,1))+Eps:L/(Nx-1):max(coords(:,1))-Eps;
B=max(coords(:,2))-min(coords(:,2))-2*Eps;
yi=min(coords(:,2))+Eps:B/(Ny-1):max(coords(:,2))-Eps;
[XI,YI]=meshgrid(xi,yi);
xp=reshape(XI,Nx*Ny,1);
yp=reshape(YI,Nx*Ny,1);
switch mode
    case 'gpu'
        ZA=zeros(NofElNodesx*NofElNodesy,Nx*Ny,'double','gpuArray');
        %ZI=zeros(Nx,Ny,'double','gpuArray');
    case 'cpu'
        ZA=zeros(NofElNodesx*NofElNodesy,Nx*Ny);
        %ZI=zeros(Nx,Ny);
end


switch surface
   case 'bottom'
   ne=1:fen;
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
bx=(x(4,:)+x(3,:))/2;
by=(y(4,:)+y(3,:))/2;
ax=(x(3,:)+x(2,:))/2;
ay=(y(3,:)+y(2,:))/2;

cx=(x(1,:)+x(4,:))/2;
cy=(y(1,:)+y(4,:))/2;
dx=(x(2,:)+x(1,:))/2;
dy=(y(2,:)+y(1,:))/2;

xc = (ax+cx)/2;
yc = (by+dy)/2;

% cross product method
% check if there is the same sign so that the point is inside convex polygon (element)
switch mode
    case 'gpu'
        ownerElement=zeros(Nx*Ny,1,'uint32','gpuArray');
        %ownerElement=zeros(Nx*Ny,1,'double','gpuArray');
    case 'cpu'
        ownerElement=zeros(Nx*Ny,1);
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

% for bounding-box test
Pmin_x=min(x);
Pmin_y=min(y);
Pmax_x=max(x);
Pmax_y=max(y);
disp('Find elements');


switch mode
    case 'gpu'
        poolobj = gcp('nocreate'); % If no pool, do not create new one.
        if isempty(poolobj)
            poolobj = parpool; % create parpool
        end
        parfor j =1:Nx*Ny
%         for j =1:Nx*Ny
            % cross product for the consequtive edge
            d1 = (xp(j).*y21-x1y21)  - (yp(j).*x21-y1x21);
            d2 = (xp(j).*y32-x2y32)  - (yp(j).*x32-y2x32);
            d3 = (xp(j).*y43-x3y43)  - (yp(j).*x43-y3x43);
            d4 = (xp(j).*y14-x4y14)  - (yp(j).*x14-y4x14);
    
            negd1=uint32(d1<0);negd2=uint32(d2<0);negd3=uint32(d3<0);negd4=uint32(d4<0);
            %negd1=(d1<0);negd2=(d2<0);negd3=(d3<0);negd4=(d4<0);
            neg=negd1.*negd2.*negd3.*negd4;
            neg_index=find(neg);
            ownerElement(j,1)=neg_index(1);
        end
        delete(poolobj);
    case 'cpu'
          % fastest option
          tic
        for j =1:Nx*Ny
            % bounding box test           
            Ibb=find(xp(j) >= Pmin_x & xp(j)<=Pmax_x & yp(j) >= Pmin_y & yp(j)<=Pmax_y);
            
            if(length(Ibb)==1)
                 ownerElement(j,1) = Ibb;
            else
                % cross product for the consequtive edges
                d1 = (xp(j).*y21(Ibb)-x1y21(Ibb))  - (yp(j).*x21(Ibb)-y1x21(Ibb));
                d2 = (xp(j).*y32(Ibb)-x2y32(Ibb))  - (yp(j).*x32(Ibb)-y2x32(Ibb));
                d3 = (xp(j).*y43(Ibb)-x3y43(Ibb))  - (yp(j).*x43(Ibb)-y3x43(Ibb));
                d4 = (xp(j).*y14(Ibb)-x4y14(Ibb))  - (yp(j).*x14(Ibb)-y4x14(Ibb));

                negd1=uint32(d1<0);negd2=uint32(d2<0);negd3=uint32(d3<0);negd4=uint32(d4<0);
                neg=negd1.*negd2.*negd3.*negd4;
                neg_index=find(neg); % can be empty
                if(isempty(neg_index))    % point is on the edge of element  
                    ownerElement(j,1)=Ibb(1);
                else
                    ownerElement(j,1)=Ibb(neg_index(1));
                end
            end
        end
        toc
end
       
disp('Calculate local values of ksi and eta for each point');
%% jacobians approach

[ksi,wx]=gll(NofElNodesx); % weights and nodes distribution
[eta,wy]=gll(NofElNodesy);

[Ksi, Eta] = meshgrid(ksi,eta);
Ksi = repmat(reshape(Ksi',1,[]),length(xp),1);
Eta = repmat(reshape(Eta',1,[]),length(xp),1);
% see: M. Li, A. Wittek, K. Miller: Efficient Inverse Isoparametric Mapping Algorithm for Whole-Body
% Computed Tomography Registration Using Deformations Predicted by Nonlinear Finite Element Modeling
% Journal of Biomechanical Engineering, Vol. 136 / 084503-1, 2014.
elementNodes_owner = nodes(ownerElement,:);
x_0 = coords(elementNodes_owner,1);
x_0 = reshape(x_0,[],NofElNodes);
y_0 = coords(elementNodes_owner,2);
y_0 = reshape(y_0,[],NofElNodes);

[~,point_no] = min(sqrt(bsxfun(@minus, x_0,xp).^2+bsxfun(@minus, y_0,yp).^2),[],2);
x0 = x_0(sub2ind(size(x_0),(1:length(point_no))',point_no));
y0 = y_0(sub2ind(size(y_0),(1:length(point_no))',point_no));

ksi0 = Ksi(sub2ind(size(Ksi),ones(length(point_no),1),point_no));
eta0 = Eta(sub2ind(size(Eta),ones(length(point_no),1),point_no));
% ksi0=repmat(ksi(1),[Nx*Ny,1]);
% eta0=repmat(eta(1),[Nx*Ny,1]);

% first iteration
n1=((ownerElement(:,1))-1)*NofElNodes+1;
ksi_p=ksi0+invj11(n1).*(xp-x0) +invj12(n1).*(yp-y0);
eta_p=eta0+invj21(n1).*(xp-x0)+invj22(n1).*(yp-y0);
clear invj11 invj12 invj21 invj22 j11 j12 j21 j22;
% n= 48451;
% ne1=ownerElement(n);
% plot(xp(n),yp(n),'ro');hold on; plot(x(:,ne1),y(:,ne1),'b');plot(x(1,ne1),y(1,ne1),'go');
[Q]=Vandermonde2D(ksi,eta,NofElNodesx,NofElNodesy);
N=shape2Dp(NofElNodesx,NofElNodesy,Q,ksi_p,eta_p);% shape functions at arbitrary (ksi, eta) point
[indxi,indxj]=find(N);
Ns=sparse(indxj,[1:NofElNodesx*NofElNodesy*Nx*Ny],reshape(N,NofElNodesx*NofElNodesy*Nx*Ny,1));
% coordinates of element nodes
xe=zeros(NofElNodesx*NofElNodesy,Nx*Ny);
ye=zeros(NofElNodesx*NofElNodesy,Nx*Ny);
for k=1:NofElNodesx*NofElNodesy
    xe(k,:)=coords(nodes((ownerElement(:,1)),k),1);
    ye(k,:)=coords(nodes((ownerElement(:,1)),k),2);
end
xe=reshape(xe,NofElNodesx*NofElNodesy*Nx*Ny,1);
ye=reshape(ye,NofElNodesx*NofElNodesy*Nx*Ny,1);
x1=Ns*xe;% interpolated values
y1=Ns*ye;% interpolated values
% plot(coords(nodes(ne1,1:NofElNodesx*NofElNodesy),1),coords(nodes(ne1,1:NofElNodesx*NofElNodesy),2),'c.');
% plot(x1(n),y1(n),'mx');
% errors
% e=sqrt((xp'-x1).^2+(yp'-y1).^2);
% [A,I]=max(e)
% J=find(e>1e-5);
% length(J)
% iterations of ksi and eta
nIterations=3;
tic
for k=1:nIterations
    ksi0=ksi_p;
    eta0=eta_p;
    x0=x1;
    y0=y1;
    % shape function derivatives
    [Npx,Npy]=shape_derivatives2Dp(NofElNodesx,NofElNodesy,Q,ksi0,eta0);
    [indxi,indxj]=find(Npx);
    Npx=sparse(indxj,[1:NofElNodesx*NofElNodesy*Nx*Ny],reshape(Npx,NofElNodesx*NofElNodesy*Nx*Ny,1));
    Npy=sparse(indxj,[1:NofElNodesx*NofElNodesy*Nx*Ny],reshape(Npy,NofElNodesx*NofElNodesy*Nx*Ny,1));
    % jacobians
    J11=Npx*xe;
    J12=Npy*xe;
    J21=Npx*ye;
    J22=Npy*ye;
    [invJ11,invJ12,invJ21,invJ22]=inv_jacp_2D(J11,J12,J21,J22);
    % next iteration
    ksi_p=ksi0+invJ11.*(xp-x0) +invJ12.*(yp-y0);
    eta_p=eta0+invJ21.*(xp-x0)+invJ22.*(yp-y0);
    
    N=shape2Dp(NofElNodesx,NofElNodesy,Q,ksi_p,eta_p);% shape functions at arbitrary (ksi, eta) point
    [indxi,indxj]=find(N);
    Ns=sparse(indxj,[1:NofElNodesx*NofElNodesy*Nx*Ny],reshape(N,NofElNodesx*NofElNodesy*Nx*Ny,1));
    x1=Ns*xe;% interpolated values
    y1=Ns*ye;% interpolated values
    %plot(x1(n),y1(n),'kx');
end
toc
% errors
e=sqrt((xp-x1).^2+(yp-y1).^2);
[A,I]=max(e);
fprintf('max error of grid points location in respect to local coordinates: %f mm \n', A*1e3);

% just in case ksi or eta are outside the range
I1=find(eta_p>1); if(I1) eta_p(I1)=1; end
I1=find(eta_p<-1); if(I1) eta_p(I1)=-1; end
I1=find(ksi_p>1); if(I1) ksi_p(I1)=1; end
I1=find(ksi_p<-1); if(I1) ksio_p(I1)=-1; end

%%
% analytic approach
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
% tic
% sx=x(1,ownerElement)+x(2,ownerElement)+x(3,ownerElement)+x(4,ownerElement);
% sy=y(1,ownerElement)+y(2,ownerElement)+y(3,ownerElement)+y(4,ownerElement);
% rx=-x(1,ownerElement)-x(2,ownerElement)+x(3,ownerElement)+x(4,ownerElement);
% ry=-y(1,ownerElement)-y(2,ownerElement)+y(3,ownerElement)+y(4,ownerElement);
% px=-x(1,ownerElement)+x(2,ownerElement)+x(3,ownerElement)-x(4,ownerElement);
% py=-y(1,ownerElement)+y(2,ownerElement)+y(3,ownerElement)-y(4,ownerElement);
% qx=x(1,ownerElement)-x(2,ownerElement)+x(3,ownerElement)-x(4,ownerElement);
% qy=y(1,ownerElement)-y(2,ownerElement)+y(3,ownerElement)-y(4,ownerElement);
% 
% A1=qy.*rx-ry.*qx;
% B1=4*yp.*qx-sy.*qx-ry.*px+py.*rx+sx.*qy-4*xp.*qy;
% C1=4*yp.*px-sy.*px-4*py.*xp+py.*sx;
% switch mode
%     case 'gpu'
%         ksi_p1 = zeros(1,Nx*Ny,'double','gpuArray');
%         eta_p1= zeros(1,Nx*Ny,'double','gpuArray');
%         ksi_p2 = zeros(1,Nx*Ny,'double','gpuArray');
%         eta_p2= zeros(1,Nx*Ny,'double','gpuArray');
%         ksi_p = zeros(1,Nx*Ny,'double','gpuArray');
%         eta_p= zeros(1,Nx*Ny,'double','gpuArray');
%         %ownerElement=zeros(Nx*Ny,1,'double','gpuArray');
%     case 'cpu'
%         ksi_p1 = zeros(1,Nx*Ny);
%         eta_p1= zeros(1,Nx*Ny);
%         ksi_p2 = zeros(1,Nx*Ny);
%         eta_p2= zeros(1,Nx*Ny);
%         ksi_p = zeros(1,Nx*Ny);
%         eta_p= zeros(1,Nx*Ny);
% end
% % needs further improvement
% II=find(A1==0); % linear equation
% eta_p1(II) = -C1(II)./B1(II);
% ksi_p1(II)=( 4*xp(II)-sx(II)-eta_p1(II).*rx(II) )./(px(II)+eta_p1(II).*qx(II));
% % resolve infinite ksi
% Iinf=find(isinf(ksi_p1));
% ksi_p1(Iinf)=( 4*yp(Iinf) -sy(Iinf)-eta_p1(Iinf).*ry(Iinf) )./( py(Iinf) + eta_p1(Iinf).*qy(Iinf));
% Iinf=find(isinf(ksi_p1));
% ksi_p1(Iinf)=0;
% % resolve nan ksi
% Inan=find(isnan(ksi_p1));
% ksi_p1(Inan) = 0;
% 
% JJ=setdiff([1:Nx*Ny],II); % quadratic equation
% Delta1=B1(JJ).^2-4*A1(JJ).*C1(JJ);
% eta_p1(JJ)= ( -B1(JJ) -sqrt(Delta1) )./(2*A1(JJ));
% eta_p2(JJ)= ( -B1(JJ) +sqrt(Delta1) )./(2*A1(JJ));
% ksi_p1(JJ)=( 4*xp(JJ)-sx(JJ)-eta_p1(JJ).*rx(JJ) )./(px(JJ)+eta_p1(JJ).*qx(JJ));
% ksi_p2(JJ)=( 4*xp(JJ)-sx(JJ)-eta_p2(JJ).*rx(JJ) )./(px(JJ)+eta_p2(JJ).*qx(JJ));
% I1=find(ksi_p1>1);
% ksi_p1(I1)=0;
% I2=find(ksi_p2>1);
% ksi_p2(I2)=0;
% I1=find(ksi_p1<-1);
% ksi_p1(I1)=0;
% I2=find(ksi_p2<-1);
% ksi_p2(I2)=0;
% % select only appropriate solutions of quadratic equation
% for k=1:Nx*Ny
%     if( abs(eta_p1(k))>1 )
%         eta_p(k) = eta_p2(k);
%         ksi_p(k) = ksi_p2(k);       
%     else
%         eta_p(k) = eta_p1(k);
%         ksi_p(k) = ksi_p1(k);
%     end
% end
% toc

% evaluate shape functions
disp('Evaluate shape functions');
% switch mode
%     case 'gpu'
%         [N]=shape_func_p_gpu(ksi_p,eta_p);
%     case 'cpu'
%         [N]=shape_func_p(ksi_p,eta_p);
% end
% % vandermonde approach
% [ksi,wx]=gll(NofElNodesx); % weights and nodes distribution
% [eta,wy]=gll(NofElNodesy);
% [Q]=Vandermonde2D(ksi,eta,NofElNodesx,NofElNodesy);
% 
% N=shape2Dp(NofElNodesx,NofElNodesy,Q,ksi_p',eta_p');
% clear ksi_p eta_p ksi_p1 eta_p1 ksi_p2 eta_p2
% % create sparse matrix
% % Ns1 = sparse(Nx*Ny,NofElNodesx*NofElNodesy*Nx*Ny);
% % for k=1:Nx*Ny
% %     Ns1(k,(k-1)*NofElNodesx*NofElNodesy+1:k*NofElNodesx*NofElNodesy) = N(:,k)';
% % end
% %
% [indxi,indxj]=find(N);
% Ns=sparse(indxj,[1:NofElNodesx*NofElNodesy*Nx*Ny],reshape(N,NofElNodesx*NofElNodesy*Nx*Ny,1));
clear N;
clear ksi_p eta_p
%%
sample=0;
for n=frm_int:frm_int:nft
    sample=sample+1;
    [sample,nft/frm_int]
    filename=fullfile('outputs',['output',num2str(k_test),],['Ux_frame',num2str(n),'.mat']);
    load(filename,'Ux'); 
    filename=fullfile('outputs',['output',num2str(k_test),],['Uy_frame',num2str(n),'.mat']);
    load(filename,'Uy'); 
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
      ZA=sqrt((Ux(nodes(ownerElement,1:NofElNodesx*NofElNodesy))').^2 + (Uy(nodes(ownerElement,1:NofElNodesx*NofElNodesy))').^2); % bottom surface
   case 'top'
      ZA=sqrt((Ux(nodes(ownerElement,NofElNodesx*NofElNodesy*nz-NofElNodesx*NofElNodesy+1:NofElNodesx*NofElNodesy*nz))').^2 + ...
          (Uy(nodes(ownerElement,NofElNodesx*NofElNodesy*nz-NofElNodesx*NofElNodesy+1:NofElNodesx*NofElNodesy*nz))').^2) ; % top surface
   case 'middle'
      ZA=sqrt((Ux(nodes(ownerElement,NofElNodesx*NofElNodesy*nz-NofElNodesx*NofElNodesy+1:NofElNodesx*NofElNodesy*nz))').^2 + ...
          (Uy(nodes(ownerElement,NofElNodesx*NofElNodesy*nz-NofElNodesx*NofElNodesy+1:NofElNodesx*NofElNodesy*nz))').^2); % middle surface
end
% c=0;
% for j =1:Ny
%     for i=1:Nx
%         c=c+1;
%         ZI(i,j) = N(:,c)'*ZA(:,c);
%     end
% end
ZA=reshape(ZA,NofElNodesx*NofElNodesy*Nx*Ny,1);
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
[i1,j1,k1] = size(Data);
if(i1*j1*k1 > 512*512*512)
    save(data_filename,'Data','-v7.3');% for data files >2GB (compression is used)
else
    save(data_filename,'Data');
end

pause(0.1);
end