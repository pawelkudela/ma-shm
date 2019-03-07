clear all;close all;   warning off;clc;
tasks=[6];
mode='gpu'; % options: mode='cpu';mode='gpu';
for k_test=tasks
    
model_output_path = fullfile('outputs',['output',num2str(k_test)],filesep);
if ~exist(model_output_path, 'dir')
    mkdir(model_output_path);
end    

assembly='mesh'; % options: assembly='trad'; assembly='mesh' ('trad' is slow);
run(fullfile('inputs',['input',num2str(k_test)]));
frm_int=floor(nft/(frames)); % save displacement with interval time step frm_int (frames)
pztEl=[];pztnum=[];

load(meshfile);

[fen,NofElNodes]=size(nodes);
[NofNodes,empty]=size(coords);

dof=5*NofNodes; % number of degrees of freedom

% host structure
[h,h1,h2,em,rhom,nim,vol,ef,rhof,nif,alpha]=lay_const(lay,lh,i_em,i_ef,i_rhom,i_rhof,i_nim,i_nif,lvol,lalpha,lmat,lfib);

clear i_em i_ef i_rhom i_rhof i_nim i_nif lvol lalpha

[ksi,wx]=gll(nx); % weights and nodes distribution
[eta,wy]=gll(ny);
[dzeta,wz]=gll(nz);
% vandermonde approach
[Qx]=Vandermonde(ksi,nx);
[Qy]=Vandermonde(eta,ny);
[Qz]=Vandermonde(dzeta,nz);
%ne=1; % element number
disp('elastic coefficient dependence on evironmental conditions')

if(isempty(glueEl)) % no glue pzt only
hoststrnum=setdiff([1:fen],pztnum);
else
pzt_bond_num=[pztnum,glueEl]; % glue and pzt
hoststrnum=setdiff([1:fen],pzt_bond_num);
end
%

%hoststrnum=1:fen; % no glue no pzt

I=zeros(fen*NofElNodes,1);
for ne=1:fen
    n1=(ne-1)*NofElNodes+1;
    n2=n1+NofElNodes-1;
    I(n1:n2,1)=nodes(ne,:);
end 
%%
% host structure
IhostG=zeros(length(hoststrnum)*NofElNodes,1);
IhostL=zeros(length(hoststrnum)*NofElNodes,1);
c=0;
for ne=hoststrnum
    c=c+1;
    n1=(c-1)*NofElNodes+1;
    n2=n1+NofElNodes-1;
    IhostG(n1:n2,1)=nodes(ne,:);
    n1L=(ne-1)*NofElNodes+1;
    n2L=n1L+NofElNodes-1;
    IhostL(n1:n2,1)=n1L:n2L;
end    


% host structure
[rho, m11, I11, a11, a22, a12, a16, a26, a66, a44, a45, a55, b11,b12,b16,b22,b26,b66,d11, d12, d16, d22, d26, d66]=composite_plate(h, h1, h2, rhom, rhof, em, ef, nim, nif, vol, alpha, lay);

%% PRELIMINARY CALCULATIONS
%disp('preliminary calculations');

%%

Inc=zeros(fen,5*NofElNodes);
% connectivity matrix Inc by dof
disp('connectivity matrix...');
 for ne=1:fen
     Inc(ne,1:5:end)=5*nodes(ne,:)-4;
     Inc(ne,2:5:end)=5*nodes(ne,:)-3;
     Inc(ne,3:5:end)=5*nodes(ne,:)-2;
     Inc(ne,4:5:end)=5*nodes(ne,:)-1;
     Inc(ne,5:5:end)=5*nodes(ne,:);
 end

%% signal
dt=tt/nft;   % calculation time step [s]
[t,st]=Hanning_signal(dt,nft,f_1,f_2,t_1);
%  figure;
%  plot(t,st);pause(1);
 
%% forces induced by the pzt actuators

%Fa2=zeros(dof,1);
%[Vx,Vzx]=Vandermonde_old(ksi,dzeta,nx,nz);
%disp('forces...');


%% Output file for solution
outfile_voltage=fullfile('outputs',['output',num2str(k_test)],['voltage',num2str(k_test)]);
outfile_displ=fullfile('outputs',['output',num2str(k_test)],['displ',num2str(k_test)]);
outfile_time=fullfile('outputs',['output',num2str(k_test)],['time',num2str(k_test)]);
%%
disp('local derivatives');
%[Nprimx,Nprimy,Nprimz]=shape3D_prim(nx,ny,nz,Qx,Qy,Qz,ksi',eta',dzeta');
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
%Z=zeros(fen*NofElNodes,1);
for ne=1:fen
    n1=(ne-1)*NofElNodes+1;
    n2=n1+NofElNodes-1;
    X(n1:n2)=coords(nodes(ne,:),1);
    Y(n1:n2)=coords(nodes(ne,:),2);
    %Z(n1:n2)=coords(nodes(ne,:),3);
end
switch mode
    case 'gpu'
    Npx=gpuArray(Npx);
    Npy=gpuArray(Npy);
 
    X=gpuArray(X);
    Y=gpuArray(Y);
    %Z=gpuArray(Z);
    case 'cpu'
        
end
disp('jacobians'); 
J11=Npx*X;
J21=Npx*Y;
%J31=Npx*Z;
J12=Npy*X;
J22=Npy*Y;
%J32=Npy*Z;


clear X Y
%clear Z
clear Indxi Indxj Indyi Indyj

%
disp('determinant and inverse Jacobian')
%%%
% /* determinant and inverse Jacobian  */

[invJ11,invJ12,invJ21,invJ22]=inv_jacp_2D(J11,J12,J21,J22);
detJ=det_jacp_2D(J11,J12,J21,J22);
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
Uc=zeros(dof,1); % displacement on cpu
switch mode
    case 'gpu'
    Iu=gpuArray(Iu);       
    U=zeros(dof,1,'double','gpuArray');
    UX=zeros(fen*NofElNodes,1,'double','gpuArray');
    UY=zeros(fen*NofElNodes,1,'double','gpuArray');
    UZ=zeros(fen*NofElNodes,1,'double','gpuArray');
    FIX=zeros(fen*NofElNodes,1,'double','gpuArray');
    FIY=zeros(fen*NofElNodes,1,'double','gpuArray');
    F2=zeros(5*fen*NofElNodes,1,'double','gpuArray');
    mgL=zeros(5*fen*NofElNodes,1,'double','gpuArray');
    mg=zeros(dof,1,'double','gpuArray'); % global mass matrix
   
    case 'cpu'
    U=zeros(dof,1);
    UX=zeros(fen*NofElNodes,1);
    UY=zeros(fen*NofElNodes,1);
    UZ=zeros(fen*NofElNodes,1); 
    FIX=zeros(fen*NofElNodes,1);
    FIY=zeros(fen*NofElNodes,1);
    F2=zeros(5*fen*NofElNodes,1);
    mgL=zeros(5*fen*NofElNodes,1);
    mg=zeros(dof,1); % global mass matrix
end

%% Pre Main       

Itemp=zeros(5*length(IG1),1);
Itemp(1:5:end)=5*IG1-4;
Itemp(2:5:end)=5*IG1-3;
Itemp(3:5:end)=5*IG1-2;
Itemp(4:5:end)=5*IG1-1;
Itemp(5:5:end)=5*IG1;
IG1=Itemp;

Itemp=zeros(5*length(IG2),1);
Itemp(1:5:end)=5*IG2-4;
Itemp(2:5:end)=5*IG2-3;
Itemp(3:5:end)=5*IG2-2;
Itemp(4:5:end)=5*IG2-1;
Itemp(5:5:end)=5*IG2;
IG2=Itemp;

Itemp=zeros(5*length(IG3),1);
Itemp(1:5:end)=5*IG3-4;
Itemp(2:5:end)=5*IG3-3;
Itemp(3:5:end)=5*IG3-2;
Itemp(4:5:end)=5*IG3-1;
Itemp(5:5:end)=5*IG3;
IG3=Itemp;

Itemp=zeros(5*length(IG4),1);
Itemp(1:5:end)=5*IG4-4;
Itemp(2:5:end)=5*IG4-3;
Itemp(3:5:end)=5*IG4-2;
Itemp(4:5:end)=5*IG4-1;
Itemp(5:5:end)=5*IG4;
IG4=Itemp;

Itemp=zeros(5*length(IG5),1);
Itemp(1:5:end)=5*IG5-4;
Itemp(2:5:end)=5*IG5-3;
Itemp(3:5:end)=5*IG5-2;
Itemp(4:5:end)=5*IG5-1;
Itemp(5:5:end)=5*IG5;
IG5=Itemp;

Itemp=zeros(5*length(IG6),1);
Itemp(1:5:end)=5*IG6-4;
Itemp(2:5:end)=5*IG6-3;
Itemp(3:5:end)=5*IG6-2;
Itemp(4:5:end)=5*IG6-1;
Itemp(5:5:end)=5*IG6;
IG6=Itemp;

Itemp=zeros(5*length(IG7),1);
Itemp(1:5:end)=5*IG7-4;
Itemp(2:5:end)=5*IG7-3;
Itemp(3:5:end)=5*IG7-2;
Itemp(4:5:end)=5*IG7-1;
Itemp(5:5:end)=5*IG7;
IG7=Itemp;

Itemp=zeros(5*length(IG8),1);
Itemp(1:5:end)=5*IG8-4;
Itemp(2:5:end)=5*IG8-3;
Itemp(3:5:end)=5*IG8-2;
Itemp(4:5:end)=5*IG8-1;
Itemp(5:5:end)=5*IG8;
IG8=Itemp;

Itemp=zeros(5*length(IG9),1);
Itemp(1:5:end)=5*IG9-4;
Itemp(2:5:end)=5*IG9-3;
Itemp(3:5:end)=5*IG9-2;
Itemp(4:5:end)=5*IG9-1;
Itemp(5:5:end)=5*IG9;
IG9=Itemp;

Itemp=zeros(5*length(IG10),1);
Itemp(1:5:end)=5*IG10-4;
Itemp(2:5:end)=5*IG10-3;
Itemp(3:5:end)=5*IG10-2;
Itemp(4:5:end)=5*IG10-1;
Itemp(5:5:end)=5*IG10;
IG10=Itemp;

Itemp=zeros(5*length(IG11),1);
Itemp(1:5:end)=5*IG11-4;
Itemp(2:5:end)=5*IG11-3;
Itemp(3:5:end)=5*IG11-2;
Itemp(4:5:end)=5*IG11-1;
Itemp(5:5:end)=5*IG11;
IG11=Itemp;

Itemp=zeros(5*length(IG12),1);
Itemp(1:5:end)=5*IG12-4;
Itemp(2:5:end)=5*IG12-3;
Itemp(3:5:end)=5*IG12-2;
Itemp(4:5:end)=5*IG12-1;
Itemp(5:5:end)=5*IG12;
IG12=Itemp;

Itemp=zeros(5*length(IL1),1);
Itemp(1:5:end)=5*IL1-4;
Itemp(2:5:end)=5*IL1-3;
Itemp(3:5:end)=5*IL1-2;
Itemp(4:5:end)=5*IL1-1;
Itemp(5:5:end)=5*IL1;
IL1=Itemp;

Itemp=zeros(5*length(IL2),1);
Itemp(1:5:end)=5*IL2-4;
Itemp(2:5:end)=5*IL2-3;
Itemp(3:5:end)=5*IL2-2;
Itemp(4:5:end)=5*IL2-1;
Itemp(5:5:end)=5*IL2;
IL2=Itemp;

Itemp=zeros(5*length(IL3),1);
Itemp(1:5:end)=5*IL3-4;
Itemp(2:5:end)=5*IL3-3;
Itemp(3:5:end)=5*IL3-2;
Itemp(4:5:end)=5*IL3-1;
Itemp(5:5:end)=5*IL3;
IL3=Itemp;

Itemp=zeros(5*length(IL4),1);
Itemp(1:5:end)=5*IL4-4;
Itemp(2:5:end)=5*IL4-3;
Itemp(3:5:end)=5*IL4-2;
Itemp(4:5:end)=5*IL4-1;
Itemp(5:5:end)=5*IL4;
IL4=Itemp;

Itemp=zeros(5*length(IL5),1);
Itemp(1:5:end)=5*IL5-4;
Itemp(2:5:end)=5*IL5-3;
Itemp(3:5:end)=5*IL5-2;
Itemp(4:5:end)=5*IL5-1;
Itemp(5:5:end)=5*IL5;
IL5=Itemp;

Itemp=zeros(5*length(IL6),1);
Itemp(1:5:end)=5*IL6-4;
Itemp(2:5:end)=5*IL6-3;
Itemp(3:5:end)=5*IL6-2;
Itemp(4:5:end)=5*IL6-1;
Itemp(5:5:end)=5*IL6;
IL6=Itemp;

Itemp=zeros(5*length(IL7),1);
Itemp(1:5:end)=5*IL7-4;
Itemp(2:5:end)=5*IL7-3;
Itemp(3:5:end)=5*IL7-2;
Itemp(4:5:end)=5*IL7-1;
Itemp(5:5:end)=5*IL7;
IL7=Itemp;

Itemp=zeros(5*length(IL8),1);
Itemp(1:5:end)=5*IL8-4;
Itemp(2:5:end)=5*IL8-3;
Itemp(3:5:end)=5*IL8-2;
Itemp(4:5:end)=5*IL8-1;
Itemp(5:5:end)=5*IL8;
IL8=Itemp;

Itemp=zeros(5*length(IL9),1);
Itemp(1:5:end)=5*IL9-4;
Itemp(2:5:end)=5*IL9-3;
Itemp(3:5:end)=5*IL9-2;
Itemp(4:5:end)=5*IL9-1;
Itemp(5:5:end)=5*IL9;
IL9=Itemp;

Itemp=zeros(5*length(IL10),1);
Itemp(1:5:end)=5*IL10-4;
Itemp(2:5:end)=5*IL10-3;
Itemp(3:5:end)=5*IL10-2;
Itemp(4:5:end)=5*IL10-1;
Itemp(5:5:end)=5*IL10;
IL10=Itemp;

Itemp=zeros(5*length(IL11),1);
Itemp(1:5:end)=5*IL11-4;
Itemp(2:5:end)=5*IL11-3;
Itemp(3:5:end)=5*IL11-2;
Itemp(4:5:end)=5*IL11-1;
Itemp(5:5:end)=5*IL11;
IL11=Itemp;

Itemp=zeros(5*length(IL12),1);
Itemp(1:5:end)=5*IL12-4;
Itemp(2:5:end)=5*IL12-3;
Itemp(3:5:end)=5*IL12-2;
Itemp(4:5:end)=5*IL12-1;
Itemp(5:5:end)=5*IL12;
IL12=Itemp;

clear Itemp;

switch mode
    case 'gpu'
    IG1=gpuArray(IG1);
    IG2=gpuArray(IG2);
    IG3=gpuArray(IG3);
    IG4=gpuArray(IG4);
    IG5=gpuArray(IG5);
    IG6=gpuArray(IG6);
    IG7=gpuArray(IG7);
    IG8=gpuArray(IG8);
    IG9=gpuArray(IG9);
    IG10=gpuArray(IG10);
    IG11=gpuArray(IG11);
    IG12=gpuArray(IG12);

    IL1=gpuArray(IL1);
    IL2=gpuArray(IL2);
    IL3=gpuArray(IL3);
    IL4=gpuArray(IL4);
    IL5=gpuArray(IL5);
    IL6=gpuArray(IL6);
    IL7=gpuArray(IL7);
    IL8=gpuArray(IL8);
    IL9=gpuArray(IL9);
    IL10=gpuArray(IL10);
    IL11=gpuArray(IL11);
    IL12=gpuArray(IL12);

    case 'cpu'
    
end
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
        mg(IG1)=mg(IG1)+mgL(IL1); 
        mg(IG2)=mg(IG2)+mgL(IL2); 
        mg(IG3)=mg(IG3)+mgL(IL3); 
        mg(IG4)=mg(IG4)+mgL(IL4); 
        mg(IG5)=mg(IG5)+mgL(IL5); 
        mg(IG6)=mg(IG6)+mgL(IL6); 
        mg(IG7)=mg(IG7)+mgL(IL7); 
        mg(IG8)=mg(IG8)+mgL(IL8); 
        mg(IG9)=mg(IG9)+mgL(IL9); 
        mg(IG10)=mg(IG10)+mgL(IL10); 
        mg(IG11)=mg(IG11)+mgL(IL11); 
        mg(IG12)=mg(IG12)+mgL(IL12); 
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
    case 'gpu'
    dv=zeros(dof,1,'double','gpuArray'); % damping vector
    case 'cpu'
    dv=zeros(dof,1); % damping vector
end

dv(1:5:end)=etad_xy; % in-plane damping factor
dv(3:5:end)=etad_xy; % in-plane damping factor
dv(5:5:end)=etad_z; % out-of-plane damping factor
cg = dv.*mg; % damping matrix
col=1;
switch mode
    case 'gpu'
    uold=zeros(dof,1,'double','gpuArray');
    unew=zeros(dof,1,'double','gpuArray'); % memory prealocation
    v=zeros(dof,1,'double','gpuArray');
    case 'cpu'
    uold=zeros(dof,1);
    unew=zeros(dof,1); % memory prealocation
    v=zeros(dof,1);
end

mg2=a2*mg;
mg=a0*mg;
mg0=1./(mg+a1*cg);

if(isempty(glueEl))
    displ=zeros(nft,length(outputs));
else
    voltage=zeros(nft,length(PZT_sensor)*length(PZT_actuator));
    displ=zeros(nft,length(outputs));
end
clear coords D V1
switch mode
    case 'gpu'
    nodes=gpuArray(nodes);
    ne=zeros(1,1,'double','gpuArray'); 
    NofElNodes=gpuArray(NofElNodes);
    n1=zeros(1,1,'double','gpuArray');
    n2=zeros(1,1,'double','gpuArray');
    case 'cpu'
    IG=[IG1;IG2;IG3;IG4;IG5;IG6;IG7;IG8;IG9;IG10;IG11;IG12];
    IL=[IL1;IL2;IL3;IL4;IL5;IL6;IL7;IL8;IL9;IL10;IL11;IL12];
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
            Fa=zeros(dof,1,'double','gpuArray');  
        case 'cpu'
            Fa=zeros(dof,1);
    end
end

if(isempty(pztEl))
    Fi(outputs(1))=1;Fi=Fi*V0; % force excitation according to dof defined in output  
    %Fi(3)=1;Fi=Fi*V0; % force excitation in z direction
else
    % pzt actuator excitation is not implemented !
end

tic;minTime = Inf;
c=0;
t_frames = zeros(frames,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('time integration...');
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
        F(IG1)=F(IG1)+F2(IL1);
        F(IG2)=F(IG2)+F2(IL2);
        F(IG3)=F(IG3)+F2(IL3);
        F(IG4)=F(IG4)+F2(IL4);
        F(IG5)=F(IG5)+F2(IL5);
        F(IG6)=F(IG6)+F2(IL6);
        F(IG7)=F(IG7)+F2(IL7);
        F(IG8)=F(IG8)+F2(IL8);
        F(IG9)=F(IG9)+F2(IL9);
        F(IG10)=F(IG10)+F2(IL10);
        F(IG11)=F(IG11)+F2(IL11);
        F(IG12)=F(IG12)+F2(IL12);
 end
    %% time integration
    %F3=-F+mg2.*U+Fa*st(nn)-mg.*uold-cg.*v;
    if(isempty(pztEl))
       % F=-F+mg2.*U+Fi*st(nn)-mg.*uold-cg.*v; % force excitation
        F=-F+mg2.*U+Fi*st(nn)-mg.*uold+a1*cg.*uold; % force excitation
    else
        %F=-F+mg2.*U+Fa*st(nn)-mg.*uold-cg.*v; % piezoelectric excitation
        F=-F+mg2.*U+Fa*st(nn)-mg.*uold+a1*cg.*uold; % piezoelectric excitation
    end  
    %
    unew=mg0.*F;
    v = a1*(-uold+unew);  % update velocity
    uold=U;  
    U = unew; 
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
    if (mod(nn,frm_int) == 0)
        c=c+1;
        t_frames(c)=t(nn);
        Uc=gather(U);
        if( max(Uc)> 10*h(1)) disp('integration error'); 
            if(isempty(pztEl))
                save(outfile_displ,'displ');
            else
                %save(outfile_voltage,'voltage');   
                save(outfile_displ,'displ'); 
            end
          return; 
        end
        if(isnan(unew(3))) disp('integration error'); 
           if(isempty(pztEl))
                save(outfile_displ,'displ');
            else
               %save(outfile_voltage,'voltage');   
               save(outfile_displ,'displ'); 
            end
           return; 
        end 
        
        outfile_Ux=fullfile('outputs',['output',num2str(k_test)],['Ux_frame' ,num2str(nn)]);
        outfile_Uy=fullfile('outputs',['output',num2str(k_test)],['Uy_frame' ,num2str(nn)]);
        outfile_Uz=fullfile('outputs',['output',num2str(k_test)],['Uz_frame' ,num2str(nn)]);
        outfile_Fix=fullfile('outputs',['output',num2str(k_test)],['Fix_frame' ,num2str(nn)]);
        outfile_Fiy=fullfile('outputs',['output',num2str(k_test)],['Fiy_frame' ,num2str(nn)]);
        
        %save(outfile,'Uc','-ascii'); % frame output  
        Ux = Uc(1:5:end); Fix = Uc(2:5:end); Uy = Uc(3:5:end); Fiy = Uc(4:5:end); Uz = Uc(5:5:end);
        save(outfile_Ux,'Ux'); % frame output for global vector of Ux displacement
        save(outfile_Uy,'Uy'); % frame output for global vector of Uy displacement
        save(outfile_Uz,'Uz'); % frame output for global vector of Uz displacement
        save(outfile_Fix,'Fix'); % frame output for global vector of Fix displacement
        save(outfile_Fiy,'Fiy'); % frame output for global vector of Fiy displacement
        averageTime = toc/(nn-1);
        progress = round(nn/nft*100);
        clc;
        message1 = sprintf('Task %d progress: %d%%',k_test,progress);
        disp(message1);
        expected_duration_seconds = (nft-nn)*averageTime; 
        s = seconds(expected_duration_seconds);
        t_now = datetime('now','TimeZone','local','Format','d-MMM-y HH:mm');
        t_expected = t_now + s;
        disp('Expected finish time: ');
        disp(t_expected);
        for jj=1:length(tasks)
            if(k_test==tasks(jj))
                break;
            end
        end
        no_of_remaining_tasks = length(tasks)-jj;
        if(no_of_remaining_tasks>0)
            remaining_tasks =tasks(jj+1:length(tasks));
            message2 = sprintf('%d remaining tasks in the queue: ',no_of_remaining_tasks);
            disp(message2);
            disp(remaining_tasks);
        else
            disp('Remaining tasks in the queue: none');
        end
    end
    telapsed = toc(tstart);
    minTime = min(telapsed,minTime);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  END OF MAIN LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
averageTime = toc/(nft-1);
%save(outfile_voltage,'voltage');   
save(outfile_displ,'displ');
t_frames_filename=fullfile('outputs',['output',num2str(k_test)],'t_frames');
save(t_frames_filename,'t_frames');
TotalNofNodes=dof/3;
save(outfile_time,'minTime','averageTime','TotalNofNodes');
switch mode
    case 'gpu'
        g = gpuDevice(1);
        reset(g);
    case 'cpu'
        
end
 
end




