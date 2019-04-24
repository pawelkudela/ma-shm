function main_flat_shell_GUI(app,AllData_filename)
% MAIN_FLAT_SHELL_GUI   Compute wavefield by using flat shell spectral elements  
% 
% Syntax: main_flat_shell_GUI(test_case,meshfile,mode,output_name)
% 
% Inputs: 
%    test_case - test case number (input/output number, variable parameters), integer
%    input_no - input file number (file with constant parameters)
%    meshfile - mesh filename, string
%    mode - string , options: mode='CPU';mode='GPU';
%    output_name - path to folder were the output data is stored, string
%    tasks - list of task numbers in a queue (for message only), integer
% 
% Outputs:
%    t_frames - time vector corresponding to frames, double, dimensions [nft,1], Units: s
% Example: 
%    main_flat_shell_GUI(test_case,meshfile,mode,output_name)
%    main_flat_shell_GUI(1,'mesh1','GPU',output_name) 
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

%load('project_paths.mat','projectroot');
test_case = 1;
tasks = [];
load(AllData_filename,'AllData');

assembly='mesh'; % options: assembly='trad'; assembly='mesh' ('trad' is slow);
model_output_path = prepare_model_paths('raw','num','flat_shell','GUI');
output_name = [model_output_path,filesep,AllData.Filename,filesep];
% check if folder exist, if not create it
if ~exist(output_name, 'dir')
    mkdir(output_name);
end
spec_mesh_output_path = fullfile(app.flat_shell_root,'mesh',filesep);
meshfile=[spec_mesh_output_path,AllData.Filename,'_mesh.mat'];
nFrames= AllData.Nooftimeframes;
nft = AllData.Noofsamples;
frm_int=floor(nft/(nFrames)); % save displacement with interval time step frm_int (frames)
%pztEl=[];pztnum=[];
mode = AllData.Mode;
load(meshfile); % coord nodes
nx=shape_order+1;
ny=shape_order+1;
[fen,NofElNodes]=size(nodes);
NofNodes=size(coords,1);
%unDelamEl=setdiff(1:fen,delamEl);
mesh_center = (mesh_max - mesh_min)/2;
[~,I] = min(sqrt( (coords(:,1) - mesh_center(1)).^2+( coords(:,2) - mesh_center(2)).^2) );
outputs=5*I; % transverse displacement at the plate center
dof=5*NofNodes; % number of degrees of freedom
app.NotificationTextArea.Value = 'Starting calculations';
drawnow;
%% selectable options
field_variable = AllData.FieldvariableDropDown;
%%
lay = AllData.lay;
lh = zeros(lay,1);
i_em = zeros(lay,1);
i_ef = zeros(lay,1);
i_rhom = zeros(lay,1);
i_rhof = zeros(lay,1);
i_nim = zeros(lay,1);
i_nif = zeros(lay,1);
lvol = zeros(lay,1);
lalpha = zeros(lay,1);
lmat = zeros(lay,1);
lfib = zeros(lay,1);

% only one material is allowed in the current implementation
etad_xy = AllData.eta_xy;
etad_z = AllData.eta_z;
% host structure material
lmat(1:16) = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];
lfib(1:16) = [1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1];
lh(1) = AllData.h1;          lh(2) = AllData.h2;         lh(3) = AllData.h3;         lh(4) = AllData.h4;
lh(5) = AllData.h5;          lh(6) = AllData.h6;         lh(7) = AllData.h7;         lh(8) = AllData.h8;
lh(9) = AllData.h9;          lh(10) = AllData.h10;      lh(11) = AllData.h11;      lh(12) = AllData.h12;
lh(13) = AllData.h13;       lh(14) = AllData.h14;      lh(15) = AllData.h15;      lh(16) = AllData.h16;
lvol(1) = AllData.vol1;     lvol(2) = AllData.vol2;     lvol(3) = AllData.vol3;     lvol(4) = AllData.vol4;
lvol(5) = AllData.vol5;     lvol(6) = AllData.vol6;     lvol(7) = AllData.vol7;     lvol(8) = AllData.vol8;
lvol(9) = AllData.vol9;     lvol(10) = AllData.vol10; lvol(11) = AllData.vol11;  lvol(12) = AllData.vol12;
lvol(13) = AllData.vol13;  lvol(14) = AllData.vol14; lvol(15) = AllData.vol15;  lvol(15) = AllData.vol15;
lalpha(1) = AllData.lay1;  lalpha(2) = AllData.lay2; lalpha(3) = AllData.lay3;  lalpha(4) = AllData.lay4;
lalpha(5) = AllData.lay5;  lalpha(6) = AllData.lay6; lalpha(7) = AllData.lay7;  lalpha(8) = AllData.lay8;
lalpha(9) = AllData.lay9;  lalpha(10) = AllData.lay10;lalpha(11) = AllData.lay11;lalpha(12) = AllData.lay12;
lalpha(13) = AllData.lay13;lalpha(14) = AllData.lay14;lalpha(15) = AllData.lay15;lalpha(16) = AllData.lay16;
i_em(1) = AllData.Em;
i_ef(1) = AllData.Ef;
i_rhom(1) = AllData.rhom;
i_rhof(1) = AllData.rhof;
i_nim(1) = AllData.nim;
i_nif(1) = AllData.nif;
lh = lh/1e3; % [m]
lvol = lvol/100;
i_em = i_em*1e9;
i_ef = i_ef*1e9;
%% pzt properties
YE11=AllData.YE11*1e9; %[N/m^2] youngs modulus at constant electric field
YE33=AllData.YE33*1e9; %[N/m^2] youngs modulus at constant electric field
nu11=AllData.nu11; % Poisson ratio

dp11=AllData.dp11;
dp12=AllData.dp12;
dp13=AllData.dp13;
dp14=AllData.dp14;
dp15=AllData.dp15;
dp16=AllData.dp16;
dp21=AllData.dp21;
dp22=AllData.dp22;
dp23=AllData.dp23;
dp24=AllData.dp24;
dp25=AllData.dp25;
dp26=AllData.dp26;
dp31=AllData.dp31;
dp32=AllData.dp32;
dp33=AllData.dp33;
dp34=AllData.dp34;
dp35=AllData.dp35;
dp36=AllData.dp36;

gp11=AllData.gp11;
gp12=AllData.gp12;
gp13=AllData.gp13;
gp14=AllData.gp14;
gp15=AllData.gp15;
gp16=AllData.gp16;
gp21=AllData.gp21;
gp22=AllData.gp22;
gp23=AllData.gp23;
gp24=AllData.gp24;
gp25=AllData.gp25;
gp26=AllData.gp26;
gp31=AllData.gp31;
gp32=AllData.gp32;
gp33=AllData.gp33;
gp34=AllData.gp34;
gp35=AllData.gp35;
gp36=AllData.gp36;
% elastic compliance matrix
 Spzt=[ 1/YE11      -nu11/YE11    -nu11/YE33        0                      0                0;
       -nu11/YE11       1/YE11       -nu11/YE33        0                      0                0;
       -nu11/YE33    -nu11/YE33     1/YE33             0                      0                0;
           0                   0                  0         (1+nu11)/YE33           0                0;
           0                   0                  0                0               (1+nu11)/YE33     0;
           0                   0                  0                0                       0         (1+nu11)/YE11];%[m^2/N] efunda
% matrix of piezoelectric coupling constants - charge constants
      %sx    sy  sz  sxz syz sxy
dp=  [ dp11 dp12 dp13 dp14 dp15 dp16;
         dp21 dp22 dp23 dp24 dp25 dp26;
         dp31 dp32 dp33 dp34 dp35 dp36]*10^-12;%[C/N] 
% voltage constants
      %sx    sy  sz  sxz syz sxy
gp=  [ gp11 gp12 gp13 gp14 gp15 gp16;
         gp21 gp22 gp23 gp24 gp25 gp26;
         gp31 gp32 gp33 gp34 gp35 gp36]*10^-3;%[Vm/N] 

% permittivity matrix (strain-charge form)
epsT=   [dp(1,5)/gp(1,5)       0                         0;
               0                 dp(1,5)/gp(1,5)            0;
               0                       0              dp(3,3)/gp(3,3)]; %[F/m]

rho_pzt=AllData.rho_pzt;%[kg/m3] % density
pzt_thickness =AllData.pzt_thickness/1e3; % pzt thickness [m]
theta_pzt = 0; % rotation angle of pzt [deg]
[h,h1,h2,em,rhom,nim,vol,ef,rhof,nif,alpha]=lay_const(lay,lh,i_em,i_ef,i_rhom,i_rhof,i_nim,i_nif,lvol,lalpha,lmat,lfib);
pause(0.01);
clear i_em i_ef i_rhom i_rhof i_nim i_nif lvol lalpha
delamination_layer = AllData.delamination_layer;
if(isempty(delamEl))
    % host structure elastic constants for the reference case
    [~, m11, I11, a11, a22, a12, a16, a26, a66, a44, a45, a55, b11,b12,b16,b22,b26,b66,d11, d12, d16, d22, d26, d66]=composite_plate(h, h1, h2, rhom, rhof, em, ef, nim, nif, vol, alpha, lay);
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
        [~, m11, I11, a11, a22, a12, a16, a26, a66, a44, a45, a55, b11,b12,b16,b22,b26,b66,d11, d12, d16, d22, d26, d66]=composite_plate_delam_pzt(h, h1, h2, rhom, rhof, em, ef, nim, nif, vol, alpha, lay, delamEl,den_above,den_under,delamination_layer,fen,NofElNodes,rho_pzt,pzt_thickness,Qpzt,pztEl);
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
tt = AllData.Totalpropagationtime/1e3;
dt=tt/nft;   % calculation time step [s]
f_1 = AllData.ModulationFrequency*1e3;% frequency of the modulation signal [Hz]
f_2 = AllData.CarrierFrequency*1e3;% frequency of the carrier signal [Hz]
t_1 = 0;
V0 =AllData.Peakvoltage; % [V]
[t,st]=Hanning_signal(dt,nft,f_1,f_2,t_1);
%   figure;
%   plot(t,st);pause(1);

%% forces induced by the pzt actuators

c=0;
PZT_actuator=[];
PZT_sensor=[];
for ne=1:1
   c=c+1;
   PZT_actuator(c).pztEl=[pztEl(ne,:)]; % list of element numbers for pzt no c
end
c=0;
for ne=1:1
   c=c+1;
   PZT_sensor(c).pztEl=[pztEl(c,:)]; % 
end



[Vx,Vzx]=Vandermonde_old(ksi,dzeta,nx,nz);

% disp('forces...');
% NofElNodes2=nx*ny*nz;
if(isempty(pztEl))
else
    Fa=zeros(dof,1);
    voltage=zeros(nft,1);
for j=1:length(PZT_actuator)
    FaxG=zeros(max(max(nodes)),1);
    FayG=zeros(max(max(nodes)),1);
    [nodespzt,coordspzt]=connect_pzt3D_to_plate_fun(nx,ny,nz,PZT_actuator(j).pztEl,nodes,coords,pzt_thickness,h);
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
    Incpzt=zeros(length(PZT_sensor(j).pztEl),3*nx*ny*nz);
    
    Incpzt(:,1:3:end)=3*nodespzt(:,:)-2;
    Incpzt(:,2:3:end)=3*nodespzt(:,:)-1;
    Incpzt(:,3:3:end)=3*nodespzt(:,:)-0;
    cpzt=0;
    for ne=PZT_actuator(j).pztEl
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
    FaL=Kfiu'*FiA; % local forces in pzt actuator j
    [boundary_nodes_pzt] = find_boundary_nodes(nodespzt(:,1:nx*ny),nx,ny);
    [boundary_nodes_global] = find_boundary_nodes(nodes(PZT_actuator(j).pztEl,1:nx*ny),nx,ny);
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
end


%% Output file for solution
outfile_voltage=fullfile(output_name,['voltage_',AllData.Filename]);
outfile_displ=fullfile(output_name,['displ_',AllData.Filename]);
outfile_time=fullfile(output_name,['time_',AllData.Filename]);
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
Xne=reshape(X,nx*ny,fen);
Yne=reshape(Y,nx*ny,fen);
switch mode
    case 'GPU'
    Npx=gpuArray(Npx);
    Npy=gpuArray(Npy);
 
    X=gpuArray(X);
    Y=gpuArray(Y);
%     Z=gpuArray(Z);
    case 'CPU'
        
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
    case 'GPU'
        j11=gather(J11);j21=gather(J21);j12=gather(J12);j22=gather(J22);
        invj11=gather(invJ11);invj21=gather(invJ21);invj12=gather(invJ12);invj22=gather(invJ22);
        det_jac=detJ;
        save([meshfile(1:end-4),'_jacobians'],'j11','j21','j12','j22','invj11','invj12','invj21','invj22','det_jac');
    case 'CPU'
        j11=J11;j21=J21;j12=J12;j22=J22;
        invj11=invJ11;invj21=invJ21;invj12=invJ12;invj22=invJ22;
        det_jac=detJ;
        save([meshfile(1:end-4),'_jacobians'],'j11','j21','j12','j22','invj11','invj12','invj21','invj22','det_jac');
end
clear j11 j12 j21 j22 invj11 invj12 invj21 invj22 det_jac
%disp('done');
%
switch mode
    case 'GPU'
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

    case 'CPU'
        
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
    case 'GPU'
    ww=gpuArray(ww);
    case 'CPU'
        
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
if(~isempty(pztEl))
    Ipzt=zeros(length(pztEl)*NofElNodes,1);
    for k=1:length(pztEl)
        ne=pztEl(k);
        n1=(ne-1)*NofElNodes+1;
        n2=n1+NofElNodes-1;
        c1=(k-1)*NofElNodes+1;
        c2=c1+NofElNodes-1;
        Ipzt(c1:c2,1)=n1:n2;
    end
end
Uc=zeros(dof,1); % displacement on CPU
switch mode
    case 'GPU'
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
   
    case 'CPU'
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
%% Initial calculation - central difference scheme
a0 = 1/(dt^2);
a1=  1/(2*dt);
a2 = 2*a0;
a3 = 1/a2;
switch mode
    case 'GPU'
    dv=zeros(dof,1,'double','gpuArray'); % damping vector
    case 'CPU'
    dv=zeros(dof,1); % damping vector
end

dv(1:5:end)=etad_xy; % in-plane damping factor
dv(3:5:end)=etad_xy; % in-plane damping factor
dv(5:5:end)=etad_z; % out-of-plane damping factor
cg = dv.*mg; % damping matrix
col=1;
switch mode
    case 'GPU'
    uold=zeros(dof,1,'double','gpuArray');
%     unew=zeros(dof,1,'double','gpuArray'); % memory prealocation
%     v=zeros(dof,1,'double','gpuArray');
    case 'CPU'
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
    case 'GPU'
    nodes=gpuArray(nodes);
    ne=zeros(1,1,'double','gpuArray'); 
    NofElNodes=gpuArray(NofElNodes);
    n1=zeros(1,1,'double','gpuArray');
    n2=zeros(1,1,'double','gpuArray');
    case 'CPU'
%     IG=[IG1;IG2;IG3;IG4;IG5;IG6;IG7;IG8;IG9;IG10;IG11;IG12];
%     IL=[IL1;IL2;IL3;IL4;IL5;IL6;IL7;IL8;IL9;IL10;IL11;IL12];
end
%
% global vector of forces induced by the actuator or nodal force excitation
if(isempty(pztEl))
    switch mode
        case 'GPU'
            Fi=zeros(dof,1,'double','gpuArray');  
        case 'CPU'
            Fi=zeros(dof,1);
    end   
else
    % pzt actuator excitation is not implemented !
    switch mode
        case 'GPU'
            Fa=gpuArray(Fa);  
        case 'CPU'
            
    end
end

if(isempty(pztEl))
    Fi(outputs(1))=1;Fi=Fi*V0; % force excitation according to dof defined in output  
    %Fi(3)=1;Fi=Fi*V0; % force excitation in z direction
end

tic;minTime = Inf;
c=0;
t_frames = zeros(nFrames,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('time integration...');
app.ActuatotrAsSensorSignal.Visible = 'on';
for nn=2:nft
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
       Volt = pzt_thickness/g33 *( e31*(NpxUX(Ipzt).*invJ11(Ipzt) +  NpyUX(Ipzt).*invJ12(Ipzt) + h/2* NpxFIX(Ipzt).*invJ11(Ipzt) + h/2*NpyFIX(Ipzt).*invJ12(Ipzt))  +...
                                             e32*(NpxUY(Ipzt).*invJ21(Ipzt) +  NpyUY(Ipzt).*invJ22(Ipzt) + h/2* NpxFIY(Ipzt).*invJ21(Ipzt) + h/2*NpyFIY(Ipzt).*invJ22(Ipzt)));
       voltage(nn,1) = gather(mean(Volt));
   end
    switch mode
        case 'GPU'
        F=zeros(dof,1,'double','gpuArray');
       
        case 'CPU'
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
                case 'GPU'
                uc1=gather(U(outputs(j1),i1));
                case 'CPU'
                uc1=U(outputs(j1),i1);
            end
            displ(nn,id)=uc1;
        end
    end
    %%
    % save frame to file
    if (mod(nn,frm_int) == 0 && c <nFrames)
        c=c+1;
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
        if( max(U)> 10*h(1)) disp('integration error'); 
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
        clc;
        
        averageTime = toc/(nn-1);
        progress = round(nn/nft*100);
        
        message1 = sprintf('Task %d progress: %d%%',test_case,progress);
        disp(message1);
        expected_duration_seconds = (nft-nn)*averageTime; 
        s = seconds(expected_duration_seconds);
        t_now = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm');
        t_expected = t_now + s;
        message2 = sprintf('Expected finish time: %s',t_expected);
        disp('Expected finish time: ');
        disp(t_expected);
        
        % plot wavefield based on corner nodes only (takes too long 25s - 30s)
%         tic
%         close all;
%         UZne = gather(UZ);
%         UZne= reshape(UZne,nx*ny,fen);
%         figure;
%         fill(Xne([1,6,36,31],:),Yne([1,6,36,31],:),UZne([1,6,36,31],:));
%         axis off ;
%         axis equal;
%         shading interp;
%         colormap jet;
%         toc
%         disp(toc);
%         
        
        
        %cla(app.ActuatotrAsSensorSignal,'reset');
        %cla(app.ActuatotrAsSensorSignal);
        
%         plot( app.ActuatotrAsSensorSignal,t(1:nn)*1e3,voltage(1:nn)*1e3 );  
%         hold(app.ActuatotrAsSensorSignal, 'on')
%         plot( app.ActuatotrAsSensorSignal,t(nn+1:end)*1e3,voltage(nn+1:end)*1e3,'r' );
        
        %plot( app.ActuatotrAsSensorSignal,t*1e3,voltage(:,1)*1e3 );
        plot( app.ActuatotrAsSensorSignal,t(1:nn)*1e3,voltage(1:nn)*1e3 ); 
        Vmax = max(abs(voltage))*1e3;
        axis(app.ActuatotrAsSensorSignal,[0 t(end)*1e3 -Vmax Vmax]);
        ylabel(app.ActuatotrAsSensorSignal,'Voltage [mV]');    
        app.NotificationTextArea.Value = message2;
        drawnow;
        if(progress >=10 && progress <20)
            app.Lamp_Progress10.Color = [0,1,0];
        end
        if(progress >=20 && progress <30)
            app.Lamp_Progress20.Color = [0,1,0];
        end
        if(progress >=30 && progress <40)
            app.Lamp_Progress30.Color = [0,1,0];
        end
        if(progress >=40 && progress <50)
            app.Lamp_Progress40.Color = [0,1,0];
        end
        if(progress >=50 && progress <60)
            app.Lamp_Progress50.Color = [0,1,0];
        end
        if(progress >=60 && progress <70)
            app.Lamp_Progress60.Color = [0,1,0];
        end
        if(progress >=70 && progress <80)
            app.Lamp_Progress70.Color = [0,1,0];
        end
        if(progress >=80 && progress <90)
            app.Lamp_Progress80.Color = [0,1,0];
        end
        if(progress >=90 && progress <100)
            app.Lamp_Progress90.Color = [0,1,0];
        end
        drawnow;
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
    end
    uold=U;  
    U = unew; 
    telapsed = toc(tstart);
    minTime = min(telapsed,minTime);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  END OF MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot( app.ActuatotrAsSensorSignal,t(1:nn)*1e3,voltage(1:nn)*1e3 ); 
Vmax = max(abs(voltage))*1e3;
axis(app.ActuatotrAsSensorSignal,[0 t(end)*1e3 -Vmax Vmax]);
ylabel(app.ActuatotrAsSensorSignal,'Voltage [mV]');    
averageTime = toc/(nft-1);
save(outfile_voltage,'voltage');   
save(outfile_displ,'displ');
model_output_path = prepare_model_paths('interim','num','flat_shell','GUI');
output_name = [model_output_path,filesep,AllData.Filename,filesep];
% check if folder exist, if not create it
if ~exist(output_name, 'dir')
    mkdir(output_name);
end
t_frames_filename=fullfile(output_name,['t_frames_',AllData.Filename]);
save(t_frames_filename,'t_frames');
TotalNofNodes=dof/3;
save(outfile_time,'minTime','averageTime','TotalNofNodes','t');
switch mode
    case 'GPU'
        g = gpuDevice(1);
        reset(g);
    case 'CPU'
        
end
t_now = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm');
message3 = sprintf('Calculation finished: %s',t_now);
app.Lamp_Progress100.Color = [0,1,0];
app.NotificationTextArea.Value = message3;
drawnow;
%---------------------- END OF CODE---------------------- 

% ================ [main_flat_shell_GUI.m] ================  
