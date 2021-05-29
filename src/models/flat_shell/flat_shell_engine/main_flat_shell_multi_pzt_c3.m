function t_frames=main_flat_shell_multi_pzt_c3(actuator_no,test_case,input_no,meshfile,mode,output_name,tasks,factor)
% main_flat_shell_multi_pzt_c3   Compute wavefield by using flat shell spectral elements and c tensor for material prop 
% 
% Syntax: main_flat_shell_multi_pzt_c3(test_case,meshfile,mode,output_name)
% 
% Inputs: 
%    test_case - test case number (input/output number, variable parameters), integer
%    input_no - input file number (file with constant parameters)
%    meshfile - mesh filename, string
%    mode - string , options: mode='cpu';mode='gpu';
%    output_name - path to folder were the output data is stored, string
%    tasks - list of task numbers in a queue (for message only), integer
%    actuator_no - pzt number corresponding to actuator, integer
%    factor - damping factor at pzt actuator location, double
% 
% Outputs:
%    t_frames - time vector corresponding to frames, double, dimensions [nft,1], Units: s
% Example: 
%    main_flat_shell_multi_pzt_c3(test_case,meshfile,mode,output_name)
%    main_flat_shell_multi_pzt_c3(1,'mesh1','gpu',output_name) 
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also: 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

assembly='mesh'; % options: assembly='trad'; assembly='mesh' ('trad' is slow);
run(fullfile('inputs',['input',num2str(input_no)]));
% save displacement with interval time step frm_int (frames)
if(mod(nft,2)) 
  % frame numbers to be saved for odd number of points as in input file
else
    frm_int = floor(nft/nFrames);
    frame_no = frm_int:frm_int:nFrames*frm_int;% frame numbers to be saved for even number of points
end

%pztEl=[];
addedMassEl = [];

load(meshfile); % coord nodes
pztnum=[];
for k=1:length(pztEl)
    pztnum=[pztnum,pztEl{k}];
end
nx=shape_order+1;
ny=shape_order+1;
[fen,NofElNodes]=size(nodes);
NofNodes=size(coords,1);
%unDelamEl=setdiff(1:fen,delamEl);
mesh_center = (mesh_max - mesh_min)/2;
[~,I] = min(sqrt( (coords(:,1) - mesh_center(1)).^2+( coords(:,2) - mesh_center(2)).^2) );
outputs=5*I; % transverse displacement at the plate center
dof=5*NofNodes; % number of degrees of freedom

% host structure material
%[h,h1,h2,em,rhom,nim,vol,ef,rhof,nif,alpha]=lay_const(lay,lh,i_em,i_ef,i_rhom,i_rhof,i_nim,i_nif,lvol,lalpha,lmat,lfib);
h2=zeros(1,lay);
h1=zeros(1,lay);
h=sum(lh(1:lay));
h2(1)=h/2; h1(1)=h/2-lh(1);
for i=2:lay
  h2(i)=h2(i-1)-lh(i-1);
  h1(i)=h1(i-1)-lh(i);
end
alpha=lalpha';
if(isempty(delamEl))
    % host structure elastic constants for the reference case
     if(isempty(pztEl))
        [~, m11, I11, a11, a22, a12, a16, a26, a66, a44, a45, a55, b11,b12,b16,b22,b26,b66,d11, d12, d16, d22, d26, d66]=composite_plate(h, h1, h2, rhom, rhof, em, ef, nim, nif, vol, alpha, lay);
     else
          % piezoelectric constants
        [Qpzt,epzt,gpzt,e31,e32,g33]=pzt_const(Spzt,dp,epsT,theta_pzt);
        % Qpzt - elastic coefficients under constant electric field
        % epzt - matrix of piezoelectric coupling constants (voltage constants)
        % gpzt - permittivity matrix in stress-charge form
        % theta_pzt - rotation angle of pzt [deg]
        if(isempty(addedMassEl))
            [~, m11, I11, a11, a22, a12, a16, a26, a66, a44, a45, a55, b11,b12,b16,b22,b26,b66,d11, d12, d16, d22, d26, d66]=composite_plate_c_delam_pzt(h, h1, h2, rho, q11,q12,q22,q44,q55,q66,alpha,lay,delamEl,den_above,den_under,delamination_layer,fen,NofElNodes,rho_pzt,pzt_thickness,Qpzt,pztnum);
            %[~, m11, I11, a11, a22, a12, a16, a26, a66, a44, a45, a55, b11,b12,b16,b22,b26,b66,d11, d12, d16, d22, d26, d66]=composite_plate_pzt(h, h1, h2, rhom, rhof, em, ef, nim, nif, vol, alpha, lay,fen,NofElNodes,rho_pzt,pzt_thickness,Qpzt,pztnum);
        else
            [~, m11, I11, a11, a22, a12, a16, a26, a66, a44, a45, a55, b11,b12,b16,b22,b26,b66,d11, d12, d16, d22, d26, d66]=composite_plate_c_added_mass_pzt(h, h1, h2, rho, q11,q12,q22,q44,q55,q66,alpha,lay,delamEl,den_above,den_under,delamination_layer,fen,NofElNodes,rho_pzt,pzt_thickness,Qpzt,pztnum,addedMassEl,rho_added_mass,added_mass_thickness,Qm);
        end
     end
else
    % host structure elastic constants for delaminated case
    if(isempty(pztEl))
        [~, m11, I11, a11, a22, a12, a16, a26, a66, a44, a45, a55, b11,b12,b16,b22,b26,b66,d11, d12, d16, d22, d26, d66]=composite_plate_delam(h, h1, h2, rhom, rhof, em, ef, nim, nif, vol, alpha, lay, delamEl,den_above,den_under,delamination_layer,fen,NofElNodes);
    else
        % piezoelectric constants
        [Qpzt,epzt,gpzt,e31,e32,g33]=pzt_const(Spzt,dp,epsT,theta_pzt);
        % Qpzt - elastic coefficients under constant electric field
        % epzt - matrix of piezoelectric coupling constants (voltage constants)
        % gpzt - permittivity matrix in stress-charge form
        % theta_pzt - rotation angle of pzt [deg]
        [~, m11, I11, a11, a22, a12, a16, a26, a66, a44, a45, a55, b11,b12,b16,b22,b26,b66,d11, d12, d16, d22, d26, d66]=composite_plate_c_delam_pzt(h, h1, h2, rho, q11,q12,q22,q44,q55,q66,alpha,lay,delamEl,den_above,den_under,delamination_layer,fen,NofElNodes,rho_pzt,pzt_thickness,Qpzt,pztnum);
      
    end
end

nz=3;

[ksi,wx]=gll(nx); % weights and nodes distribution
[eta,wy]=gll(ny);
[dzeta,wz]=gll(nz);
% vandermonde approach
[Qx]=Vandermonde(ksi,nx);
[Qy]=Vandermonde(eta,ny);
[Qz]=Vandermonde(dzeta,nz);

I=zeros(fen*NofElNodes,1);
for ne=1:fen
    n1=(ne-1)*NofElNodes+1;
    n2=n1+NofElNodes-1;
    I(n1:n2,1)=nodes(ne,:);
end 
%%
% host structure
IhostG=zeros(fen*NofElNodes,1);
IhostL=zeros(fen*NofElNodes,1);
c=0;
for ne=1:fen
    c=c+1;
    n1=(c-1)*NofElNodes+1;
    n2=n1+NofElNodes-1;
    IhostG(n1:n2,1)=nodes(ne,:);
    n1L=(ne-1)*NofElNodes+1;
    n2L=n1L+NofElNodes-1;
    IhostL(n1:n2,1)=n1L:n2L;
end    

%% PRELIMINARY CALCULATIONS
%disp('preliminary calculations');

%% signal
dt=tt/(nft-1);   % calculation time step [s]
[t,st]=Hanning_signal(dt,nft,f_1,f_2,t_1);
%  figure;
%  plot(t,st);pause(1);
 
%% forces induced by the pzt actuators

PZT_actuator = pztEl{actuator_no}; % only one transducer act as actuator
PZT_sensor = pztEl; % all transducers are sensors
[Vx,Vzx]=Vandermonde_old(ksi,dzeta,nx,nz);

% disp('forces...');
% NofElNodes2=nx*ny*nz;
if(isempty(pztEl))
else
    Fa=zeros(dof,1);
    voltage=zeros(nft,length(PZT_sensor));

    FaxG=zeros(max(max(nodes)),1);
    FayG=zeros(max(max(nodes)),1);
    [nodespzt,coordspzt]=connect_pzt3D_to_plate_fun(nx,ny,nz,PZT_actuator,nodes,coords,pzt_thickness,h);
    %     [p,nodespzt]=remove_free_spec_nodes(coords,nodespzt);
    % electric boundary conditions for sensor in open circuit
    
    EBC0=nodespzt(:,1:nx*ny);
    EBCI=nodespzt(:,nx*ny+1:end-nx*ny);
    EBCN=nodespzt(:,end-nx*ny+1:end);
    
    EBC0=unique(EBC0);
    EBCI=unique(EBCI);
    EBCN=unique(EBCN);
    EBC0=reshape(EBC0,1,[]);
    EBCI=reshape(EBCI,1,[]);
    EBCN=reshape(EBCN,1,[]);
    BC03=zeros(1,3*length(EBC0));
    BC03(1:3:end)=3*EBC0-2;
    BC03(2:3:end)=3*EBC0-1;
    BC03(3:3:end)=3*EBC0-0;
    
    BCI3=zeros(1,3*length(EBCI));
    BCI3(1:3:end)=3*EBCI-2;
    BCI3(2:3:end)=3*EBCI-1;
    BCI3(3:3:end)=3*EBCI-0;
    
    BCN3=zeros(1,3*length(EBCN));
    BCN3(1:3:end)=3*EBCN-2;
    BCN3(2:3:end)=3*EBCN-1;
    BCN3(3:3:end)=3*EBCN-0;
    % assembly
    NofNodes = max(max(nodespzt));
    dofpzt=NofNodes*3;
    Kfifi=sparse(NofNodes,NofNodes);
    Kfiu=sparse(NofNodes,dofpzt);
    Incpzt=zeros(length(PZT_actuator),3*nx*ny*nz);
    
    Incpzt(:,1:3:end)=3*nodespzt(:,:)-2;
    Incpzt(:,2:3:end)=3*nodespzt(:,:)-1;
    Incpzt(:,3:3:end)=3*nodespzt(:,:)-0;
    cpzt=0;
    for ne=PZT_actuator
        cpzt=cpzt+1;
        [kufi,kfifi]=coupl3D_global(epzt,gpzt,coordspzt(nodespzt(cpzt,:),1),coordspzt(nodespzt(cpzt,:),2),(coordspzt(nodespzt(cpzt,:),3)),Qx,Qy,Qz,ksi,eta,dzeta,wx,wy,wz);
        
        kfiu=sparse(kufi');kfifi=sparse(kfifi);
        [iout]=dofs3D(Incpzt,cpzt,dofpzt,nx*ny*nz);
        [iout2]=dofs3Dfifi(nodespzt,cpzt,NofNodes,nx*ny*nz);
        Kfiu=Kfiu+iout2'*kfiu*iout;
        Kfifi=Kfifi+iout2'*kfifi*iout2;
    end
    KfiuA=[Kfiu(EBCI,BC03),Kfiu(EBCI,BCI3),Kfiu(EBCI,BCN3)];% actuator BC
    KfifiA=Kfifi(EBCI,EBCI); % actuator
    Kfifiin=Kfifi(EBCI,EBCN);
    V1=sparse(NofNodes,1);
    V1(EBCN,1)=V0; % top
    V1(EBC0,1)=0; % bottom
    Fetemp=Kfifiin*V1(EBCN);
    fiAi=-KfifiA\Fetemp; % induced potential
    FiA=sparse(NofNodes,1);
    
    FiA(EBCN)=V0;
    FiA(EBCI)=fiAi;
    FaL=Kfiu'*FiA; % local forces in pzt actuator
    [boundary_nodes_pzt] = find_boundary_nodes(nodespzt(:,1:nx*ny),nx,ny);
    [boundary_nodes_global] = find_boundary_nodes(nodes(PZT_actuator,1:nx*ny),nx,ny);
    FaxL=FaL(1:3:end);
    FayL=FaL(2:3:end);
    FazL=FaL(3:3:end);
    FaxG(boundary_nodes_global,1)=FaxL(boundary_nodes_pzt,1);
    FayG(boundary_nodes_global,1)=FayL(boundary_nodes_pzt,1);
    Fa(1:5:end,1) = FaxG; % pzt equivalent force along x axis
    Fa(3:5:end,1) = FayG; % pzt equivalent force alogn y axis
    Fa(2:5:end,1) = FaxG *(h/2+pzt_thickness/2); % pzt equivalent bending moment at neutral plane
    Fa(4:5:end,1) = FayG *(h/2+pzt_thickness/2); % pzt equivalent bending moment at neutral plane

end

scale = 1;

%% Output file for solution
outfile_voltage=fullfile(output_name,['voltage',num2str(test_case)]);
outfile_displ=fullfile(output_name,['displ',num2str(test_case)]);
outfile_time=fullfile(output_name,['time',num2str(test_case)]);
%%
disp('local derivatives');
[Nprimx,Nprimy]=shape2D_prim(nx,ny,Qx,Qy,ksi',eta');
%[N]=shape2D(nx,ny,Qx,Qy,ksi',eta');

[indxi,indxj]=find(Nprimx);% indices of nonzero values for derivative in x direction
[indyi,indyj]=find(Nprimy);% indices of nonzero values for derivative in y direction

ixl=length(indxi);
iyl=length(indyi);

Npxs=zeros(ixl,1);
Npys=zeros(iyl,1);

for k=1:ixl
    Npxs(k)=Nprimx(indxi(k),indxj(k));      
end
for k=1:iyl
    Npys(k)=Nprimy(indyi(k),indyj(k));
end

Nxs=repmat(Npxs,fen,1);
Nys=repmat(Npys,fen,1);

clear Npxs Npys Npzs 
Indxi=zeros(fen*ixl,1);
Indxj=zeros(fen*ixl,1);
Indyi=zeros(fen*iyl,1);
Indyj=zeros(fen*iyl,1);

for ne=1:fen
    offset=(ne-1)*NofElNodes;
    n1=(ne-1)*ixl+1;
    n2=n1+ixl-1;
    n3=(ne-1)*iyl+1;
    n4=n3+iyl-1;
    Indxi(n1:n2,1)=offset+indxi;
    Indxj(n1:n2,1)=offset+indxj;
    Indyi(n3:n4,1)=offset+indyi;
    Indyj(n3:n4,1)=offset+indyj;
end
Npx=sparse(Indxi,Indxj,Nxs);
clear Nxs;
Npy=sparse(Indyi,Indyj,Nys);
clear Nys;

clear indxi indxj indyi indyj;
X=zeros(fen*NofElNodes,1);
Y=zeros(fen*NofElNodes,1);
% Z=zeros(fen*NofElNodes,1);
for ne=1:fen
    n1=(ne-1)*NofElNodes+1;
    n2=n1+NofElNodes-1;
    X(n1:n2)=coords(nodes(ne,:),1);
    Y(n1:n2)=coords(nodes(ne,:),2);
%     Z(n1:n2)=coords(nodes(ne,:),3);
end
switch mode
    case 'gpu'
    Npx=gpuArray(Npx);
    Npy=gpuArray(Npy);
 
    X=gpuArray(X);
    Y=gpuArray(Y);
%     Z=gpuArray(Z);
    case 'cpu'
        
end
disp('jacobians'); 
J11=Npx*X;
J21=Npx*Y;
% J31=Npx*Z;
J12=Npy*X;
J22=Npy*Y;
% J32=Npy*Z;



% clear X Y
% clear Z
clear Indxi Indxj Indyi Indyj

%
disp('determinant and inverse Jacobian')
%%%
% /* determinant and inverse Jacobian  */

[invJ11,invJ12,invJ21,invJ22]=inv_jacp_2D(J11,J12,J21,J22);
detJ=det_jacp_2D(J11,J12,J21,J22);
switch mode
    case 'gpu'
        j11=gather(J11);j21=gather(J21);j12=gather(J12);j22=gather(J22);
        invj11=gather(invJ11);invj21=gather(invJ21);invj12=gather(invJ12);invj22=gather(invJ22);
        det_jac=gather(detJ);
        save([meshfile,'_jacobians'],'j11','j21','j12','j22','invj11','invj12','invj21','invj22','det_jac');
    case 'cpu'
        j11=J11;j21=J21;j12=J12;j22=J22;
        invj11=invJ11;invj21=invJ21;invj12=invJ12;invj22=invJ22;
        det_jac=detJ;
        save([meshfile,'_jacobians'],'j11','j21','j12','j22','invj11','invj12','invj21','invj22','det_jac');
end
fprintf('Min determinant of Jacobi matrix: %e\n',min(det_jac));
fprintf('Max determinant of Jacobi matrix: %e\n',max(det_jac));
clear j11 j12 j21 j22 invj11 invj12 invj21 invj22 det_jac
%disp('done');
%
switch mode
    case 'gpu'
        a44=gpuArray(a44);
        a45=gpuArray(a45);
        a55=gpuArray(a55);
  
        a11=gpuArray(a11);
        a12=gpuArray(a12);
        a16=gpuArray(a16);
        a22=gpuArray(a22);
        a26=gpuArray(a26);
        a66=gpuArray(a66);
  
        b11=gpuArray(b11);
        b12=gpuArray(b12);
        b16=gpuArray(b16);
        b22=gpuArray(b22);
        b26=gpuArray(b26);
        b66=gpuArray(b66);
  
        d11=gpuArray(d11);
        d12=gpuArray(d12);
        d16=gpuArray(d16);
        d22=gpuArray(d22);
        d26=gpuArray(d26);
        d66=gpuArray(d66);

    case 'cpu'
        
end
%%
      
% weights at quadrature points
cc=0;
ww=zeros(NofElNodes,1);

    for ky=1:ny % eta
        for kx=1:nx % ksi  
            cc=cc+1;
            ww(cc,1)=wx(kx)*wy(ky);
        end
    end

switch mode
    case 'gpu'
    ww=gpuArray(ww);
    case 'cpu'
        
end
% 
WW=repmat(ww,fen,1);
WWDetJ=WW.*detJ;

clear WW ww detJ
clear J11 J12 J21 J22 J31 J32
  
Iu=zeros(fen*NofElNodes,1);
Ifix=zeros(fen*NofElNodes,1);
Iv=zeros(fen*NofElNodes,1);
Ifiy=zeros(fen*NofElNodes,1);
Iw=zeros(fen*NofElNodes,1);
for ne=1:fen
    n1=(ne-1)*NofElNodes+1;
    n2=n1+NofElNodes-1;
    Iu(n1:n2,1)=5*nodes(ne,:)-4;
    Ifix(n1:n2,1)=5*nodes(ne,:)-3;
    Iv(n1:n2,1)=5*nodes(ne,:)-2;
    Ifiy(n1:n2,1)=5*nodes(ne,:)-1;
    Iw(n1:n2,1)=5*nodes(ne,:);
end    
% degrees of freedom to acquire sensor signals
if(~isempty(pztEl))
    Ipzt=cell(1,length(PZT_sensor));
    for j =1:length(PZT_sensor)
        sensor_elements = PZT_sensor{j};
        Isensor = zeros(1,length(sensor_elements)*NofElNodes);
        for k=1:length(sensor_elements)
            ne=sensor_elements(k);
            n1=(ne-1)*NofElNodes+1;
            n2=n1+NofElNodes-1;
            c1=(k-1)*NofElNodes+1;
            c2=c1+NofElNodes-1;
            Isensor(c1:c2)=n1:n2;
        end
        Ipzt{j} = Isensor;
    end
end
Uc=zeros(dof,1); % displacement on cpu
switch mode
    case 'gpu'
    Iu=gpuArray(Iu);       
    U=zeros(dof,1,'double','gpuArray');
%     UX=zeros(fen*NofElNodes,1,'double','gpuArray');
%     UY=zeros(fen*NofElNodes,1,'double','gpuArray');
%     UZ=zeros(fen*NofElNodes,1,'double','gpuArray');
%     FIX=zeros(fen*NofElNodes,1,'double','gpuArray');
%     FIY=zeros(fen*NofElNodes,1,'double','gpuArray');
    F2=zeros(5*fen*NofElNodes,1,'double','gpuArray');
    mgL=zeros(5*fen*NofElNodes,1,'double','gpuArray');
    mg=zeros(dof,1,'double','gpuArray'); % global mass matrix
   
    case 'cpu'
    U=zeros(dof,1);
%     UX=zeros(fen*NofElNodes,1);
%     UY=zeros(fen*NofElNodes,1);
%     UZ=zeros(fen*NofElNodes,1); 
%     FIX=zeros(fen*NofElNodes,1);
%     FIY=zeros(fen*NofElNodes,1);
    F2=zeros(5*fen*NofElNodes,1);
    mgL=zeros(5*fen*NofElNodes,1);
    mg=zeros(dof,1); % global mass matrix
end

%% Pre Main       
%[IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12] = global_local_12baskets_5dofs(IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12,mode);
[I_G_dofs,I_L_dofs] = global_local_12baskets_5dofs_v2(I_G,I_L,mode);
%% mass matrix
disp('assembly mass matrix');
%

switch assembly
     case 'trad'
        for ne=1:fen
             n1=(ne-1)*NofElNodes+1;
             n2=n1+NofElNodes-1;              
             mg(5*nodes(ne,:)-4,1)=mg(5*nodes(ne,:)-4,1)+WWDetJ(n1:n2).*m11;
             mg(5*nodes(ne,:)-3,1)=mg(5*nodes(ne,:)-3,1)+WWDetJ(n1:n2).*I11;
             mg(5*nodes(ne,:)-2,1)=mg(5*nodes(ne,:)-2,1)+WWDetJ(n1:n2).*m11;
             mg(5*nodes(ne,:)-1,1)=mg(5*nodes(ne,:)-1,1)+WWDetJ(n1:n2).*I11;
             mg(5*nodes(ne,:),1)=mg(5*nodes(ne,:),1)+WWDetJ(n1:n2).*m11;
        end 
    case 'mesh'   
        mgL(1:5:end)=WWDetJ.*m11; % mass matrix (global node numbering -x direction)
        mgL(2:5:end)=WWDetJ.*I11; % mass matrix (global node numbering -fix direction)
        mgL(3:5:end)=WWDetJ.*m11; % mass matrix (global node numbering -y direction)
        mgL(4:5:end)=WWDetJ.*I11; % mass matrix (global node numbering -fiy direction)
        mgL(5:5:end)=WWDetJ.*m11; % mass matrix (global node numbering -z direction)
        for ibasket=1:12
            mg(I_G_dofs(:,ibasket))=mg(I_G_dofs(:,ibasket))+mgL(I_L_dofs(:,ibasket)); 
        end
end 
clear m mgL;
%
%
%% Boundary conditions ??
%  apply to Npx, Npy, ?? Npxt
%%
Npxt=Npx';
Npyt=Npy';
%%

%% Initial calculation - central difference scheme
a0 = 1/(dt^2);
a1=  1/(2*dt);
a2 = 2*a0;
a3 = 1/a2;
switch mode
    case 'gpu'
    dv=zeros(dof,1,'double','gpuArray'); % damping vector
    case 'cpu'
    dv=zeros(dof,1); % damping vector
end

dv(1:5:end)=etad_xy; % in-plane damping factor
dv(2:5:end)=etad_z2; % out-of-plane damping factor; rotational dof
dv(3:5:end)=etad_xy; % in-plane damping factor
dv(4:5:end)=etad_z2; % out-of-plane damping factor; rotational dof
dv(5:5:end)=etad_z; % out-of-plane damping factor

% higher damping at actuator location
pzt_nodes=unique(reshape(nodes(PZT_actuator,:),[],1));
dv(5*pzt_nodes-4)=factor*etad_xy; % in-plane damping factor
dv(5*pzt_nodes-3)=factor*etad_z2; % out-of-plane damping factor; rotational dof
dv(5*pzt_nodes-2)=factor*etad_xy; % in-plane damping factor
dv(5*pzt_nodes-1)=factor*etad_z2; % out-of-plane damping factor; rotational dof
dv(5*pzt_nodes-0)=factor*etad_z; % out-of-plane damping factor
cg = dv.*mg; % damping matrix
col=1;
switch mode
    case 'gpu'
    uold=zeros(dof,1,'double','gpuArray');
%     unew=zeros(dof,1,'double','gpuArray'); % memory prealocation
%     v=zeros(dof,1,'double','gpuArray');
    case 'cpu'
    uold=zeros(dof,1);
%     unew=zeros(dof,1); % memory prealocation
%     v=zeros(dof,1);
end

mg2=a2*mg;
mg=a0*mg;
mg0=1./(mg+a1*cg);

displ=zeros(nft,length(outputs));

clear coords D V1
switch mode
    case 'gpu'
    nodes=gpuArray(nodes);
    ne=zeros(1,1,'double','gpuArray'); 
    NofElNodes=gpuArray(NofElNodes);
    n1=zeros(1,1,'double','gpuArray');
    n2=zeros(1,1,'double','gpuArray');
    case 'cpu'
%     IG=[IG1;IG2;IG3;IG4;IG5;IG6;IG7;IG8;IG9;IG10;IG11;IG12];
%     IL=[IL1;IL2;IL3;IL4;IL5;IL6;IL7;IL8;IL9;IL10;IL11;IL12];
end
%
% global vector of forces induced by the actuator or nodal force excitation
if(isempty(pztEl))
    switch mode
        case 'gpu'
            Fi=zeros(dof,1,'double','gpuArray');  
        case 'cpu'
            Fi=zeros(dof,1);
    end   
else
    % pzt actuator excitation is not implemented !
    switch mode
        case 'gpu'
            Fa=gpuArray(Fa);  
        case 'cpu'
            
    end
end

if(isempty(pztEl))
    Fi(outputs(1))=1;Fi=Fi*V0; % force excitation according to dof defined in output  
    %Fi(3)=1;Fi=Fi*V0; % force excitation in z direction
end

tic;minTime = Inf;
c=1;
t_frames = zeros(nFrames,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('time integration...');
for nn=1:nft
    tstart = tic;
    %[nn,nft]
% transform from global displacement vector to local elemental vector
    UX=U(Iu);  
    UY=U(Iv);
    UZ=U(Iw);
    FIX=U(Ifix);
    FIY=U(Ifiy);
    
    NpxUX=Npx*UX;
    NpxUY=Npx*UY;
    NpxFIX=Npx*FIX;
    NpxFIY=Npx*FIY;
    
    NpyUX=Npy*UX;
    NpyUY=Npy*UY;
    NpyFIX=Npy*FIX;
    NpyFIY=Npy*FIY;
    NpxUZ=Npx*UZ;
    NpyUZ=Npy*UZ;
    
    % only valid for flat plate (z=0) isoparametric element
   Sxx= ((NpxUX).*invJ11+(NpyUX).*invJ21).*a11+((NpxFIX).*invJ11+(NpyFIX).*invJ21).*b11+((NpxUY).*invJ12+(NpyUY).*invJ22).*a12+((NpxFIY).*invJ12+(NpyFIY).*invJ22).*b12+((NpxUX).*invJ12+(NpyUX).*invJ22).*a16+((NpxFIX).*invJ12+(NpyFIX).*invJ22).*b16+((NpxUY).*invJ11+(NpyUY).*invJ21).*a16+((NpxFIY).*invJ11+(NpyFIY).*invJ21).*b16;
   Sxxz=((NpxUX).*invJ11+(NpyUX).*invJ21).*b11+((NpxFIX).*invJ11+(NpyFIX).*invJ21).*d11+((NpxUY).*invJ12+(NpyUY).*invJ22).*b12+((NpxFIY).*invJ12+(NpyFIY).*invJ22).*d12+((NpxUX).*invJ12+(NpyUX).*invJ22).*b16+((NpxFIX).*invJ12+(NpyFIX).*invJ22).*d16+((NpxUY).*invJ11+(NpyUY).*invJ21).*b16+((NpxFIY).*invJ11+(NpyFIY).*invJ21).*d16;
   Sxy= ((NpxUX).*invJ11+(NpyUX).*invJ21).*a16+((NpxFIX).*invJ11+(NpyFIX).*invJ21).*b16+((NpxUY).*invJ12+(NpyUY).*invJ22).*a26+((NpxFIY).*invJ12+(NpyFIY).*invJ22).*b26+((NpxUX).*invJ12+(NpyUX).*invJ22).*a66+((NpxFIX).*invJ12+(NpyFIX).*invJ22).*b66+((NpxUY).*invJ11+(NpyUY).*invJ21).*a66+((NpxFIY).*invJ11+(NpyFIY).*invJ21).*b66;
   Sxyz=((NpxUX).*invJ11+(NpyUX).*invJ21).*b16+((NpxFIX).*invJ11+(NpyFIX).*invJ21).*d16+((NpxUY).*invJ12+(NpyUY).*invJ22).*b26+((NpxFIY).*invJ12+(NpyFIY).*invJ22).*d26+((NpxUX).*invJ12+(NpyUX).*invJ22).*b66+((NpxFIX).*invJ12+(NpyFIX).*invJ22).*d66+((NpxUY).*invJ11+(NpyUY).*invJ21).*b66+((NpxFIY).*invJ11+(NpyFIY).*invJ21).*d66;
   Syy= ((NpxUX).*invJ11+(NpyUX).*invJ21).*a12+((NpxFIX).*invJ11+(NpyFIX).*invJ21).*b12+((NpxUY).*invJ12+(NpyUY).*invJ22).*a22+((NpxFIY).*invJ12+(NpyFIY).*invJ22).*b22+((NpxUX).*invJ12+(NpyUX).*invJ22).*a26+((NpxFIX).*invJ12+(NpyFIX).*invJ22).*b26+((NpxUY).*invJ11+(NpyUY).*invJ21).*a26+((NpxFIY).*invJ11+(NpyFIY).*invJ21).*b26;
   Syyz=((NpxUX).*invJ11+(NpyUX).*invJ21).*b12+((NpxFIX).*invJ11+(NpyFIX).*invJ21).*d12+((NpxUY).*invJ12+(NpyUY).*invJ22).*b22+((NpxFIY).*invJ12+(NpyFIY).*invJ22).*d22+((NpxUX).*invJ12+(NpyUX).*invJ22).*b26+((NpxFIX).*invJ12+(NpyFIX).*invJ22).*d26+((NpxUY).*invJ11+(NpyUY).*invJ21).*b26+((NpxFIY).*invJ11+(NpyFIY).*invJ21).*d26;
   Syz=(FIY+(NpxUZ).*invJ12+(NpyUZ).*invJ22).*a44+(FIX+(NpxUZ).*invJ11+(NpyUZ).*invJ21).*a45;
   Sxz=(FIY+(NpxUZ).*invJ12+(NpyUZ).*invJ22).*a45+(FIX+(NpxUZ).*invJ11+(NpyUZ).*invJ21).*a55;
   
   fu=Npxt*(Sxx.*invJ11.*WWDetJ)+Npyt*(Sxx.*invJ21.*WWDetJ)+Npxt*(Sxy.*invJ12.*WWDetJ)+Npyt*(Sxy.*invJ22.*WWDetJ);  
   fv=Npxt*(Syy.*invJ12.*WWDetJ)+Npyt*(Syy.*invJ22.*WWDetJ)+Npxt*(Sxy.*invJ11.*WWDetJ)+Npyt*(Sxy.*invJ21.*WWDetJ);

   mx=Npxt*(Sxxz.*invJ11.*WWDetJ)+Npyt*(Sxxz.*invJ21.*WWDetJ)+Npxt*(Sxyz.*invJ12.*WWDetJ)+Npyt*(Sxyz.*invJ22.*WWDetJ)+Sxz.*WWDetJ;
   my=Npxt*(Syyz.*invJ12.*WWDetJ)+Npyt*(Syyz.*invJ22.*WWDetJ)+Npxt*(Sxyz.*invJ11.*WWDetJ)+Npyt*(Sxyz.*invJ21.*WWDetJ)+Syz.*WWDetJ;
   
   fw=Npxt*(Syz.*invJ12.*WWDetJ)+Npyt*(Syz.*invJ22.*WWDetJ)+Npxt*(Sxz.*invJ11.*WWDetJ)+Npyt*(Sxz.*invJ21.*WWDetJ);
   if(~isempty(pztEl))
       for k=1:length(PZT_sensor)
           Volt = pzt_thickness/g33 *( e31*(NpxUX(Ipzt{k}).*invJ11(Ipzt{k}) +  NpyUX(Ipzt{k}).*invJ12(Ipzt{k}) + h/2* NpxFIX(Ipzt{k}).*invJ11(Ipzt{k}) + h/2*NpyFIX(Ipzt{k}).*invJ12(Ipzt{k}))  +...
                                                 e32*(NpxUY(Ipzt{k}).*invJ21(Ipzt{k}) +  NpyUY(Ipzt{k}).*invJ22(Ipzt{k}) + h/2* NpxFIY(Ipzt{k}).*invJ21(Ipzt{k}) + h/2*NpyFIY(Ipzt{k}).*invJ22(Ipzt{k})));
           voltage(nn,k) = gather(mean(Volt));
       end
   end
    switch mode
        case 'gpu'
        F=zeros(dof,1,'double','gpuArray');
       
        case 'cpu'
        F=zeros(dof,1);
    end
    
    
 switch assembly
     case 'trad'
        for ne=1:fen
             n1=(ne-1)*NofElNodes+1;
             n2=n1+NofElNodes-1;              
             F(5*nodes(ne,:)-4,1)=F(5*nodes(ne,:)-4,1)+fu(n1:n2);
             F(5*nodes(ne,:)-3,1)=F(5*nodes(ne,:)-3,1)+mx(n1:n2);
             F(5*nodes(ne,:)-2,1)=F(5*nodes(ne,:)-2,1)+fv(n1:n2);
             F(5*nodes(ne,:)-1,1)=F(5*nodes(ne,:)-1,1)+my(n1:n2);
             F(5*nodes(ne,:),1)=F(5*nodes(ne,:),1)+fw(n1:n2);
        end 
        clear fu fv fw   
     case 'mesh'  
        F2(1:5:end)=fu;
        F2(2:5:end)=mx;
        F2(3:5:end)=fv;
        F2(4:5:end)=my;
        F2(5:5:end)=fw;
        clear fu fv fw
        %
        for ibasket=1:12
            F(I_G_dofs(:,ibasket))=F(I_G_dofs(:,ibasket))+F2(I_L_dofs(:,ibasket)); 
        end
%         F(IG1)=F(IG1)+F2(IL1);
%         F(IG2)=F(IG2)+F2(IL2);
%         F(IG3)=F(IG3)+F2(IL3);
%         F(IG4)=F(IG4)+F2(IL4);
%         F(IG5)=F(IG5)+F2(IL5);
%         F(IG6)=F(IG6)+F2(IL6);
%         F(IG7)=F(IG7)+F2(IL7);
%         F(IG8)=F(IG8)+F2(IL8);
%         F(IG9)=F(IG9)+F2(IL9);
%         F(IG10)=F(IG10)+F2(IL10);
%         F(IG11)=F(IG11)+F2(IL11);
%         F(IG12)=F(IG12)+F2(IL12);
 end
    %% time integration
    %F3=-F+mg2.*U+Fa*st(nn)-mg.*uold-cg.*v;
    if(isempty(pztEl))
        %F=-F+mg2.*U+Fi*st(nn)-mg.*uold-cg.*v; % force excitation
        F=-F+mg2.*U+Fi*st(nn)-mg.*uold+a1*cg.*uold; % force excitation
    else
        %F=-F+mg2.*U+Fa*st(nn)-mg.*uold-cg.*v; % piezoelectric excitation
        F=-F+mg2.*U+Fa*st(nn)-mg.*uold+a1*cg.*uold; % piezoelectric excitation
    end  
    %
    unew=mg0.*F;
    %% output signals
   id=0;
    for i1=1:col
        for j1=1:length(outputs)
        id=id+1;
            switch mode
                case 'gpu'
                uc1=gather(U(outputs(j1),i1));
                case 'cpu'
                uc1=U(outputs(j1),i1);
            end
            displ(nn,id)=uc1;
        end
    end
    %%
    % save frame to file
    if (nn==frame_no(c)) 
        
        t_frames(c)=t(nn);    
        switch field_variable
            case 'displacement'
                Uc=gather(U);
            case 'velocity'
                v = a1*(-uold+unew);  % update velocity
                Vc=gather(v);
            case 'acceleration'
                anew = a0*(uold-2*U+unew); % update acceleration
                Ac=gather(anew);
            case 'all'
                Uc=gather(U);
                v = a1*(-uold+unew);  % update velocity
                Vc=gather(v);
                anew = a0*(uold-2*U+unew); % update acceleration
                Ac=gather(anew);
        end
        if( max(UZ)> 10*h(1)) disp('integration error'); 
            disp(max(UZ))
            if(isempty(pztEl))
                save(outfile_displ,'displ');
            else
                save(outfile_voltage,'voltage');   
                save(outfile_displ,'displ'); 
            end
          return; 
        end
        if(isnan(unew(3))) disp('integration error'); 
           if(isempty(pztEl))
                save(outfile_displ,'displ');
            else
               save(outfile_voltage,'voltage');   
               save(outfile_displ,'displ'); 
            end
           return; 
        end 
        % %    field_variable - string: 'displacement', 'velocity', 'acceleration' or 'all'
        switch field_variable
            case 'displacement'
                outfile_Ux=fullfile(output_name,['Ux_frame' ,num2str(c)]);
                outfile_Uy=fullfile(output_name,['Uy_frame' ,num2str(c)]);
                outfile_Uz=fullfile(output_name,['Uz_frame' ,num2str(c)]);
                outfile_UFix=fullfile(output_name,['UFix_frame' ,num2str(c)]);
                outfile_UFiy=fullfile(output_name,['UFiy_frame' ,num2str(c)]);

                Ux = Uc(1:5:end); UFix = Uc(2:5:end); Uy = Uc(3:5:end); UFiy = Uc(4:5:end); Uz = Uc(5:5:end);
                save(outfile_Ux,'Ux'); % frame output for global vector of Ux displacement
                save(outfile_Uy,'Uy'); % frame output for global vector of Uy displacement
                save(outfile_Uz,'Uz'); % frame output for global vector of Uz displacement
                save(outfile_UFix,'UFix'); % frame output for global vector of Fix displacement
                save(outfile_UFiy,'UFiy'); % frame output for global vector of Fiy displacement
            case 'velocity'      
                outfile_Vx=fullfile(output_name,['Vx_frame' ,num2str(c)]);
                outfile_Vy=fullfile(output_name,['Vy_frame' ,num2str(c)]);
                outfile_Vz=fullfile(output_name,['Vz_frame' ,num2str(c)]);
                outfile_VFix=fullfile(output_name,['VFix_frame' ,num2str(c)]);
                outfile_VFiy=fullfile(output_name,['VFiy_frame' ,num2str(c)]);
                
                Vx = Vc(1:5:end); VFix = Vc(2:5:end); Vy = Vc(3:5:end); VFiy = Vc(4:5:end); Vz = Vc(5:5:end);
                save(outfile_Vx,'Vx'); % frame output for global vector of Vx velocity
                save(outfile_Vy,'Vy'); % frame output for global vector of Vy velocity
                save(outfile_Vz,'Vz'); % frame output for global vector of Vz velocity
                save(outfile_VFix,'VFix'); % frame output for global vector of VFix velocity
                save(outfile_VFiy,'VFiy'); % frame output for global vector of VFiy velocity
            case 'acceleration'
                outfile_Ax=fullfile(output_name,['Ax_frame' ,num2str(c)]);
                outfile_Ay=fullfile(output_name,['Ay_frame' ,num2str(c)]);
                outfile_Az=fullfile(output_name,['Az_frame' ,num2str(c)]);
                outfile_AFix=fullfile(output_name,['AFix_frame' ,num2str(c)]);
                outfile_AFiy=fullfile(output_name,['AFiy_frame' ,num2str(c)]);

                Ax = Ac(1:5:end); AFix = Ac(2:5:end); Ay = Ac(3:5:end); AFiy = Ac(4:5:end); Az = Ac(5:5:end);
                save(outfile_Ax,'Ax'); % frame output for global vector of Ax acceleration
                save(outfile_Ay,'Ay'); % frame output for global vector of Ay acceleration
                save(outfile_Az,'Az'); % frame output for global vector of Az acceleration
                save(outfile_AFix,'AFix'); % frame output for global vector of AFix acceleration
                save(outfile_AFiy,'AFiy'); % frame output for global vector of AFiy acceleration
            case 'all'
                outfile_Ux=fullfile(output_name,['Ux_frame' ,num2str(c)]);
                outfile_Uy=fullfile(output_name,['Uy_frame' ,num2str(c)]);
                outfile_Uz=fullfile(output_name,['Uz_frame' ,num2str(c)]);
                outfile_UFix=fullfile(output_name,['UFix_frame' ,num2str(c)]);
                outfile_UFiy=fullfile(output_name,['UFiy_frame' ,num2str(c)]);

                Ux = Uc(1:5:end); UFix = Uc(2:5:end); Uy = Uc(3:5:end); UFiy = Uc(4:5:end); Uz = Uc(5:5:end);
                save(outfile_Ux,'Ux'); % frame output for global vector of Ux displacement
                save(outfile_Uy,'Uy'); % frame output for global vector of Uy displacement
                save(outfile_Uz,'Uz'); % frame output for global vector of Uz displacement
                save(outfile_UFix,'UFix'); % frame output for global vector of Fix displacement
                save(outfile_UFiy,'UFiy'); % frame output for global vector of Fiy displacement
                
                outfile_Vx=fullfile(output_name,['Vx_frame' ,num2str(c)]);
                outfile_Vy=fullfile(output_name,['Vy_frame' ,num2str(c)]);
                outfile_Vz=fullfile(output_name,['Vz_frame' ,num2str(c)]);
                outfile_VFix=fullfile(output_name,['VFix_frame' ,num2str(c)]);
                outfile_VFiy=fullfile(output_name,['VFiy_frame' ,num2str(c)]);
                
                Vx = Vc(1:5:end); VFix = Vc(2:5:end); Vy = Vc(3:5:end); VFiy = Vc(4:5:end); Vz = Vc(5:5:end);
                save(outfile_Vx,'Vx'); % frame output for global vector of Vx velocity
                save(outfile_Vy,'Vy'); % frame output for global vector of Vy velocity
                save(outfile_Vz,'Vz'); % frame output for global vector of Vz velocity
                save(outfile_VFix,'VFix'); % frame output for global vector of VFix velocity
                save(outfile_VFiy,'VFiy'); % frame output for global vector of VFiy velocity
                
                outfile_Ax=fullfile(output_name,['Ax_frame' ,num2str(c)]);
                outfile_Ay=fullfile(output_name,['Ay_frame' ,num2str(c)]);
                outfile_Az=fullfile(output_name,['Az_frame' ,num2str(c)]);
                outfile_AFix=fullfile(output_name,['AFix_frame' ,num2str(c)]);
                outfile_AFiy=fullfile(output_name,['AFiy_frame' ,num2str(c)]);

                Ax = Ac(1:5:end); AFix = Ac(2:5:end); Ay = Ac(3:5:end); AFiy = Ac(4:5:end); Az = Ac(5:5:end);
                save(outfile_Ax,'Ax'); % frame output for global vector of Ax acceleration
                save(outfile_Ay,'Ay'); % frame output for global vector of Ay acceleration
                save(outfile_Az,'Az'); % frame output for global vector of Az acceleration
                save(outfile_AFix,'AFix'); % frame output for global vector of AFix acceleration
                save(outfile_AFiy,'AFiy'); % frame output for global vector of AFiy acceleration   
        end
        averageTime = toc/(nn-1);
        progress = round(nn/nft*100);
        clc;
        message1 = sprintf('Task %d progress: %d%%',test_case,progress);
        disp(message1);
        expected_duration_seconds = (nft-nn)*averageTime; 
        s = seconds(expected_duration_seconds);
        t_now = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm');
        t_expected = t_now + s;
        disp('Expected finish time: ');
        disp(t_expected);
        for jj=1:length(tasks)
            if(test_case==tasks(jj))
                break;
            end
        end
        no_of_remaining_tasks = length(tasks)-jj;
        if(no_of_remaining_tasks>0)
            %remaining_tasks =tasks(jj+1:length(tasks));
             if(no_of_remaining_tasks==1)
                message2 = sprintf('%d remaining task in the queue ',no_of_remaining_tasks);
             else
                 message2 = sprintf('%d remaining tasks in the queue ',no_of_remaining_tasks);
             end
            disp(message2);
            %disp(remaining_tasks);
        else
            disp('Remaining tasks in the queue: none');
        end
        if(c<length(frame_no)) c=c+1; end
    end
    uold=U;  
    U = unew; 
    telapsed = toc(tstart);
    minTime = min(telapsed,minTime);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  END OF MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
averageTime = toc/(nft-1);
save(outfile_voltage,'voltage');   
save(outfile_displ,'displ');
% t_frames_filename=fullfile(output_name,'t_frames');
% save(t_frames_filename,'t_frames');
TotalNofNodes=dof/3;
save(outfile_time,'minTime','averageTime','TotalNofNodes','t');
switch mode
    case 'gpu'
        g = gpuDevice(1);
        reset(g);
    case 'cpu'
        
end

%---------------------- END OF CODE---------------------- 

% ================ [main_flat_shell_multi_pzt_c.m] ================  
