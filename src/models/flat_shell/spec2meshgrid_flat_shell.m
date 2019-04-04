function [Data] = spec2meshgrid_flat_shell(test_case,meshfile,Nx,Ny,field_variable,motion,shell_surface,input_name,output_name)
% SPEC2MESHGRID_FLAT_SHELL   Interpolate wavefield on uniform grid 
%    wavefield is spanned on a spectral element mesh 
%    uniform mesh is created by using meshgrid  
% 
% Syntax: [Data] = spec2meshgrid_flat_shell(test_case,meshfile,Nx,Ny,field_variable,motion,shell_surface,output_name) 
% 
% Inputs: 
%    test_case - test case number (input/output number), integer
%    meshfile - mesh filename, string
%    Nx - Number of points in uniform mesh in x direction, integer 
%    Ny - Number of points in uniform mesh in y direction, integer 
%    field_variable - string: 'displacement', 'velocity' or 'acceleration' 
%    motion - integer defining motion type: 
%    1) Ux
%    2) Uy
%    3) Uz
%    4) Fix
%    5) Fiy
%    6) Ux+h/2*UFix
%    7) Uy+h/2*UFiy
%    8) sqrt((Ux+h/2.*UFix).^2+(Uy+h/2.*UFiy).^2)
%    9) sqrt((Ux+h/2.*UFix).^2+(Uy+h/2.*UFiy).^2 + Uz.^2)
%    shell_surface - string: 'top' or 'bottom', field variable interpolated by using nodes above or under 
%    delamination and +h/2, -h/2, respectively
%    input_name - path to folder were the input data is stored, string
%    output_name - path to folder were the output data is stored, string
%
% Outputs: 
%    Data - Interpolated wavefield, double, dimensions [Nx, Ny, numberOfTimeFrames], Units: m, m/s or  m/s^2 
% 
% Example: 
%    [Data] = spec2meshgrid_flat_shell(test_case,meshfile,Nx,Ny,field_variable,motion,shell_surface) 
%    [Data] = spec2meshgrid_flat_shell(test_case,meshfile,Nx,Ny,'displacement',1,'top') 
%    [Data] = spec2meshgrid_flat_shell(test_case,meshfile,Nx,Ny,'velocity',4,'bottom') 
% 
% Other m-files required: gll.m flat_shell_variable_names.m shape2Dp.m shape_derivatives2Dp.m
%                                Vandermonde2D.m
% Subfunctions: none 
% MAT-files required:  Uz_frame*.mat Ux_frame*.mat Uy_frame*.mat Fix_frame*.mat Fiy_frame*.mat
%                            Vz_frame*.mat Vx_frame*.mat Vy_frame*.mat VFix_frame*.mat VFiy_frame*.mat
%                            Az_frame*.mat Ax_frame*.mat Ay_frame*.mat AFix_frame*.mat AFiy_frame*.mat
% See also: 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

k_test=test_case;
Eps=1e-8;
input_file = fullfile('inputs',['input',num2str(k_test)]);
run(input_file); 
load([meshfile(1:end-4),'_jacobians']);
load(meshfile); % cords, nodes
h=sum(lh);

%delamnum=[];
[variable_name] = flat_shell_variable_names(field_variable,motion);

%data_filename=fullfile('outputs',['\output',num2str(k_test)],['\flat_shell_',variable_name,'_',num2str(k_test),'_',num2str(Nx),'x',num2str(Ny),shell_surface,'.mat']);
data_filename=fullfile(output_name,['\flat_shell_',variable_name,'_',num2str(k_test),'_',num2str(Nx),'x',num2str(Ny),shell_surface,'.mat']);
% corner nodes
NofElNodesx = shape_order +1;
NofElNodesy = shape_order +1;
mc=[1,NofElNodesx,NofElNodesx*NofElNodesy,NofElNodesx*NofElNodesy-NofElNodesx+1];
% % corner nodes
% switch NofElNodesx
%      case 10
%         mc=[1,10,100,91];
%     case 9
%         mc=[1,9,81,73];
%      case 8
%         mc=[1,8,64,57];
%     case 7
%         mc=[1,7,49,43];
%     case 6 % 36-nodes element
%         mc=[1,6,36,31];
%     case 5
%         mc=[1,5,25,21];
%     case 4
%         mc=[1,4,16,13];
%     case 3
%         mc=[1,3,9,7];
% end
Data=zeros(Ny,Nx,nFrames);
[fen,NofElNodes]=size(nodes);
unDelamEl=setdiff(1:fen,delamEl);
if(isempty(delamEl))
    ne=1:fen;
else
    switch shell_surface
        case 'top'
            for k=1:length(den_above)
                fen=fen-length(den_above{k});
            end
        case 'bottom'
            for k=1:length(den_under)
                fen=fen-length(den_under{k});
            end 
    end
    ne=1:fen;
end

L=max(coords(:,1))-min(coords(:,1))-2*Eps;
xi=min(coords(:,1))+Eps:L/(Nx-1):max(coords(:,1))-Eps;
B=max(coords(:,2))-min(coords(:,2))-2*Eps;
yi=min(coords(:,2))+Eps:B/(Ny-1):max(coords(:,2))-Eps;
[XI,YI]=meshgrid(xi,yi);
xp=reshape(XI,1,Nx*Ny);
yp=reshape(YI,1,Nx*Ny);
%ZA=zeros(NofElNodesx*NofElNodesy,Nx*Ny);



%%
% extract corner nodes of all elements
x=coords(nodes(ne,mc),1);
y=coords(nodes(ne,mc),2);


x=reshape(x,fen,4)';
y=reshape(y,fen,4)';
 
% % find centres (xc,yc) of all elements
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

elem_index=zeros(Nx*Ny,1);

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
disp('Determination of the owner element');

for j =1:Nx*Ny
    % bounding box test           
    Ibb=find(xp(j) >= Pmin_x & xp(j)<=Pmax_x & yp(j) >= Pmin_y & yp(j)<=Pmax_y);

    if(length(Ibb)==1)
         elem_index(j,1) = Ibb;
    else
        % cross product for the consequtive edges
        d1 = (xp(j).*y21(Ibb)-x1y21(Ibb))  - (yp(j).*x21(Ibb)-y1x21(Ibb));
        d2 = (xp(j).*y32(Ibb)-x2y32(Ibb))  - (yp(j).*x32(Ibb)-y2x32(Ibb));
        d3 = (xp(j).*y43(Ibb)-x3y43(Ibb))  - (yp(j).*x43(Ibb)-y3x43(Ibb));
        d4 = (xp(j).*y14(Ibb)-x4y14(Ibb))  - (yp(j).*x14(Ibb)-y4x14(Ibb));

        negd1=uint32(d1<0);negd2=uint32(d2<0);negd3=uint32(d3<0);negd4=uint32(d4<0);
        neg=negd1.*negd2.*negd3.*negd4;
        neg_index=find(neg);
        elem_index(j,1)=Ibb(neg_index(1));
    end
end
        
       
disp('Calculate local values of ksi and eta for each point');
%% jacobians approach

[ksi,wx]=gll(NofElNodesx); % weights and nodes distribution
[eta,wy]=gll(NofElNodesy);

% see: M. Li, A. Wittek, K. Miller: Efficient Inverse Isoparametric Mapping Algorithm for Whole-Body
% Computed Tomography Registration Using Deformations Predicted by Nonlinear Finite Element Modeling
% Journal of Biomechanical Engineering, Vol. 136 / 084503-1, 2014.
ksi0=repmat(ksi(1),[Nx*Ny,1]);
eta0=repmat(eta(1),[Nx*Ny,1]);
x0=x(1,elem_index(:,1))';
y0=y(1,elem_index(:,1))';
% first iteration
n1=((elem_index(:,1))-1)*NofElNodes+1;
ksi_p=ksi0+invj11(n1).*(xp'-x0) +invj12(n1).*(yp'-y0);
eta_p=eta0+invj21(n1).*(xp'-x0)+invj22(n1).*(yp'-y0);
clear invj11 invj12 invj21 invj22 j11 j12 j21 j22;
% n= 48451;
% ne1=elem_index(n);
% plot(xp(n),yp(n),'ro');hold on; plot(x(:,ne1),y(:,ne1),'b');plot(x(1,ne1),y(1,ne1),'go');
[Q]=Vandermonde2D(ksi,eta,NofElNodesx,NofElNodesy);
N=shape2Dp(NofElNodesx,NofElNodesy,Q,ksi_p,eta_p);% shape functions at arbitrary (ksi, eta) point
[indxi,indxj]=find(N);
Ns=sparse(indxj,[1:NofElNodesx*NofElNodesy*Nx*Ny],reshape(N,NofElNodesx*NofElNodesy*Nx*Ny,1));
% coordinates of element nodes
xe=zeros(NofElNodesx*NofElNodesy,Nx*Ny);
ye=zeros(NofElNodesx*NofElNodesy,Nx*Ny);
for k=1:NofElNodesx*NofElNodesy
    xe(k,:)=coords(nodes((elem_index(:,1)),k),1);
    ye(k,:)=coords(nodes((elem_index(:,1)),k),2);
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
    ksi_p=ksi0+invJ11.*(xp'-x0) +invJ12.*(yp'-y0);
    eta_p=eta0+invJ21.*(xp'-x0)+invJ22.*(yp'-y0);
    
    N=shape2Dp(NofElNodesx,NofElNodesy,Q,ksi_p,eta_p);% shape functions at arbitrary (ksi, eta) point
    [indxi,indxj]=find(N);
    Ns=sparse(indxj,[1:NofElNodesx*NofElNodesy*Nx*Ny],reshape(N,NofElNodesx*NofElNodesy*Nx*Ny,1));
    x1=Ns*xe;% interpolated values
    y1=Ns*ye;% interpolated values
    %plot(x1(n),y1(n),'kx');
end

% errors
e=sqrt((xp'-x1).^2+(yp'-y1).^2);
[A,I]=max(e);
fprintf('max error of grid points location in respect to local coordinates: %f mm \n', A*1e3);

% just in case ksi or eta are outside the range
I1=find(eta_p>1); if(I1) eta_p(I1)=1; end
I1=find(eta_p<-1); if(I1) eta_p(I1)=-1; end
I1=find(ksi_p>1); if(I1) ksi_p(I1)=1; end
I1=find(ksi_p<-1); if(I1) ksio_p(I1)=-1; end

% evaluate shape functions
disp('Evaluate shape functions');
clear N;
clear ksi_p eta_p
switch shell_surface
    case 'top'
        nodes_above = [];
        nodes_under = [];
        for k=1:length(den_above)
            nodes_above = [nodes_above,nodes(den_above{k},:)];     
        end
        for k=1:length(den_under)
            nodes_under = [nodes_under,nodes(den_under{k},:)];     
        end
end
%%
for n=1:nFrames
    switch field_variable
        case 'displacement'
            switch motion
                case 1
                    filename=fullfile(input_name,['Ux_frame',num2str(n),'.mat']);
                    load(filename,'Ux'); 
                    U=Ux;
                case 2
                    filename=fullfile(input_name,['Uy_frame',num2str(n),'.mat']);
                    load(filename,'Uy'); 
                    U=Uy;
                case 3
                    filename=fullfile(input_name,['Uz_frame',num2str(n),'.mat']);
                    load(filename,'Uz'); 
                    U=Uz;
                case 4
                    filename=fullfile(input_name,['UFix_frame',num2str(n),'.mat']);
                    load(filename,'UFix'); 
                    U=UFix;
                case 5
                    filename=fullfile(input_name,['UFiy_frame',num2str(n),'.mat']);
                    load(filename,'UFiy'); 
                    U=UFiy;
                case 6
                    filename=fullfile(input_name,['Ux_frame',num2str(n),'.mat']);
                    load(filename,'Ux'); 
                    filename=fullfile(input_name,['UFix_frame',num2str(n),'.mat']);
                    load(filename,'UFix'); 
                    switch shell_surface
                        case 'top'
                            U=Ux+h/2.*UFix;
                        case 'bottom'
                            U=Ux-h/2.*UFix;
                    end
                case 7
                    filename=fullfile(input_name,['Uy_frame',num2str(n),'.mat']);
                    load(filename,'Uy'); 
                    filename=fullfile(input_name,['UFiy_frame',num2str(n),'.mat']);
                    load(filename,'UFiy'); 
                    switch shell_surface
                        case 'top'
                            U=Uy+h/2.*UFiy;
                        case 'bottom'
                            U=Uy-h/2.*UFiy;
                    end
                case 8
                    filename=fullfile(input_name,['Ux_frame',num2str(n),'.mat']);
                    load(filename,'Ux'); 
                    filename=fullfile(input_name,['UFix_frame',num2str(n),'.mat']);
                    load(filename,'UFix'); 
                    filename=fullfile(input_name,['Uy_frame',num2str(n),'.mat']);
                    load(filename,'Uy'); 
                    filename=fullfile(input_name,['UFiy_frame',num2str(n),'.mat']);
                    load(filename,'UFiy'); 
                    switch shell_surface
                        case 'top'
                            U=sqrt((Ux+h/2.*UFix).^2+(Uy+h/2.*UFiy).^2);
                        case 'bottom'
                            U=sqrt((Ux-h/2.*UFix).^2+(Uy-h/2.*UFiy).^2);
                    end
                case 9
                    filename=fullfile(input_name,['Ux_frame',num2str(n),'.mat']);
                    load(filename,'Ux'); 
                    filename=fullfile(input_name,['UFix_frame',num2str(n),'.mat']);
                    load(filename,'UFix'); 
                    filename=fullfile(input_name,['Uy_frame',num2str(n),'.mat']);
                    load(filename,'Uy'); 
                    filename=fullfile(input_name,['UFiy_frame',num2str(n),'.mat']);
                    load(filename,'UFiy'); 
                    filename=fullfile(input_name,['Uz_frame',num2str(n),'.mat']);
                    load(filename,'Uz'); 
                    switch shell_surface
                        case 'top'
                            U=sqrt((Ux+h/2.*UFix).^2+(Uy+h/2.*UFiy).^2 + Uz.^2);
                        case 'bottom'
                            U=sqrt((Ux-h/2.*UFix).^2+(Uy-h/2.*UFiy).^2 + Uz.^2);
                    end
            end
        case 'velocity'
            switch motion
                case 1
                    filename=fullfile(input_name,['Vx_frame',num2str(n),'.mat']);
                    load(filename,'Vx'); 
                    U=Vx;
                case 2
                    filename=fullfile(input_name,['Vy_frame',num2str(n),'.mat']);
                    load(filename,'Vy'); 
                    U=Vy;
                case 3
                    filename=fullfile(input_name,['Vz_frame',num2str(n),'.mat']);
                    load(filename,'Vz'); 
                    U=Vz;
                case 4
                    filename=fullfile(input_name,['VFix_frame',num2str(n),'.mat']);
                    load(filename,'VFix'); 
                    U=VFix;
                case 5
                    filename=fullfile(input_name,['VFiy_frame',num2str(n),'.mat']);
                    load(filename,'VFiy'); 
                    U=VFiy;
                case 6
                    filename=fullfile(input_name,['Vx_frame',num2str(n),'.mat']);
                    load(filename,'Vx'); 
                    filename=fullfile(input_name,['VFix_frame',num2str(n),'.mat']);
                    load(filename,'VFix'); 
                    switch shell_surface
                        case 'top'
                            U=Vx+h/2.*VFix;
                        case 'bottom'
                             U=Vx-h/2.*VFix;
                    end
                case 7
                    filename=fullfile(input_name,['Vy_frame',num2str(n),'.mat']);
                    load(filename,'Vy'); 
                    filename=fullfile(input_name,['VFiy_frame',num2str(n),'.mat']);
                    load(filename,'VFiy'); 
                    switch shell_surface
                        case 'top'
                            U=Vy+h/2.*VFiy;
                        case 'bottom'
                            U=Vy-h/2.*VFiy;
                    end
                case 8
                    filename=fullfile(input_name,['Vx_frame',num2str(n),'.mat']);
                    load(filename,'Vx'); 
                    filename=fullfile(input_name,['VFix_frame',num2str(n),'.mat']);
                    load(filename,'VFix'); 
                    filename=fullfile(input_name,['Vy_frame',num2str(n),'.mat']);
                    load(filename,'Vy'); 
                    filename=fullfile(input_name,['VFiy_frame',num2str(n),'.mat']);
                    load(filename,'VFiy'); 
                    switch shell_surface
                        case 'top'
                            U=sqrt((Vx+h/2.*VFix).^2+(Vy+h/2.*VFiy).^2);
                        case 'bottom'
                            U=sqrt((Vx-h/2.*VFix).^2+(Vy-h/2.*VFiy).^2);
                    end
                case 9
                    filename=fullfile(input_name,['Vx_frame',num2str(n),'.mat']);
                    load(filename,'Vx'); 
                    filename=fullfile(input_name,['VFix_frame',num2str(n),'.mat']);
                    load(filename,'VFix'); 
                    filename=fullfile(input_name,['Vy_frame',num2str(n),'.mat']);
                    load(filename,'Vy'); 
                    filename=fullfile(input_name,['VFiy_frame',num2str(n),'.mat']);
                    load(filename,'VFiy'); 
                    filename=fullfile(input_name,['Vz_frame',num2str(n),'.mat']);
                    load(filename,'Vz'); 
                    switch shell_surface
                        case 'top'
                            U=sqrt((Vx+h/2.*VFix).^2+(Vy+h/2.*VFiy).^2 + Vz.^2);
                        case 'bottom'
                            U=sqrt((Vx-h/2.*VFix).^2+(Vy-h/2.*VFiy).^2 + Vz.^2);
                    end
            end
        case 'acceleration'
            switch motion
                case 1
                    filename=fullfile(input_name,['Ax_frame',num2str(n),'.mat']);
                    load(filename,'Ax'); 
                    U=Ax;
                case 2
                    filename=fullfile(input_name,['Ay_frame',num2str(n),'.mat']);
                    load(filename,'Ay'); 
                    U=Ay;
                case 3
                    filename=fullfile(input_name,['Az_frame',num2str(n),'.mat']);
                    load(filename,'Az'); 
                    U=Az;
                case 4
                    filename=fullfile(input_name,['AFix_frame',num2str(n),'.mat']);
                    load(filename,'AFix'); 
                    U=AFix;
                case 5
                    filename=fullfile(input_name,['AFiy_frame',num2str(n),'.mat']);
                    load(filename,'AFiy'); 
                    U=AFiy;
                case 6
                    filename=fullfile(input_name,['Ax_frame',num2str(n),'.mat']);
                    load(filename,'Ax'); 
                    filename=fullfile(input_name,['AFix_frame',num2str(n),'.mat']);
                    load(filename,'AFix'); 
                    switch shell_surface
                        case 'top'
                            U=Ax+h/2.*AFix;
                        case 'bottom'
                            U=Ax-h/2.*AFix;
                    end
                case 7
                    filename=fullfile(input_name,['Ay_frame',num2str(n),'.mat']);
                    load(filename,'Ay'); 
                    filename=fullfile(input_name,['AFiy_frame',num2str(n),'.mat']);
                    load(filename,'AFiy'); 
                    switch shell_surface
                        case 'top'
                            U=Ay+h/2.*AFiy;
                        case 'bottom'
                            U=Ay-h/2.*AFiy;
                    end
                case 8
                    filename=fullfile(input_name,['Ax_frame',num2str(n),'.mat']);
                    load(filename,'Ax'); 
                    filename=fullfile(input_name,['AFix_frame',num2str(n),'.mat']);
                    load(filename,'AFix'); 
                    filename=fullfile(input_name,['Ay_frame',num2str(n),'.mat']);
                    load(filename,'Ay'); 
                    filename=fullfile(input_name,['AFiy_frame',num2str(n),'.mat']);
                    load(filename,'AFiy'); 
                    switch shell_surface
                        case 'top'
                            U=sqrt((Ax+h/2.*AFix).^2+(Ay+h/2.*AFiy).^2);
                        case 'bottom'
                            U=sqrt((Ax-h/2.*AFix).^2+(Ay-h/2.*AFiy).^2);
                    end
                case 9
                    filename=fullfile(input_name,['Ax_frame',num2str(n),'.mat']);
                    load(filename,'Ax'); 
                    filename=fullfile(input_name,['AFix_frame',num2str(n),'.mat']);
                    load(filename,'AFix'); 
                    filename=fullfile(input_name,['Ay_frame',num2str(n),'.mat']);
                    load(filename,'Ay'); 
                    filename=fullfile(input_name,['AFiy_frame',num2str(n),'.mat']);
                    load(filename,'AFiy'); 
                    filename=fullfile(input_name,['Az_frame',num2str(n),'.mat']);
                    load(filename,'Az'); 
                    switch shell_surface
                        case 'top'
                            U=sqrt((Ax+h/2.*AFix).^2+(Ay+h/2.*AFiy).^2 + Az.^2);
                        case 'bottom'
                            U=sqrt((Ax-h/2.*AFix).^2+(Ay-h/2.*AFiy).^2 + Az.^2);
                    end
            end
    end
    switch shell_surface
        case 'top'
            U(nodes_under) = U(nodes_above); % overwrite values of bottom part of delamination by upper part
    end
    ZA=U(nodes(elem_index,1:NofElNodesx*NofElNodesy))'; 
    ZA=reshape(ZA,NofElNodesx*NofElNodesy*Nx*Ny,1);
    ZI=Ns*ZA;% interpolated values
    ZI=reshape(ZI,Nx,Ny);
    Data(:,:,n)=ZI;

end

save(data_filename,'Data');

%---------------------- END OF CODE---------------------- 

% ================ [spec2meshgrid_flat_shell.m] ================  
