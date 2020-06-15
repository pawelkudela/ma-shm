clear all;close all;clc;warning off; clc;
tasks = [7:8];
mode='gpu'; % options: mode='cpu';mode='gpu';
for k_test=tasks

%% make sure that those lines are executed in bash
% export PATH="$PATH:~/usr/local/jacket/engine"
% export LD_LIBRARY_PATH=/usr/local/jacket/engine:$LD_LIBRARY_PATH
model_output_path = fullfile('outputs',['output',num2str(k_test)],filesep);
if ~exist(model_output_path, 'dir')
    mkdir(model_output_path);
end    

assembly='mesh'; % options: assembly='trad'; assembly='mesh' ('trad' is slow);
multi='multipzt'; % options: multi='multifocus';multi='multipzt';
run(fullfile('inputs',['input',num2str(k_test)]));

% Preliminary calculations
eval(['load ', meshfile]); % load mesh: two matrices coords and nodes and vectors IL1:IL12 IG1:IG12
coords = coords3D;
nodes = nodes3D;
[fen,NofElNodes]=size(nodes);
[NofNodes,empty]=size(coords);
%eval(['load ', envfile]); % load environmental effects: two vectors Temp and Moist
dof=3*NofNodes; % number of degrees of freedom
Temp(1:dof,1)=20; Moist(1:dof,1)=1; % dry
% piezoelectric constants
[Qpzt,epzt,gpzt]=pzt_const(Spzt,dp,epsT);
% Qpzt - elastic coefficients under constant electric field
% epzt - matrix of piezoelectric coupling constants (voltage constants)
% gpzt - permittivity matrix in stress-charge form

% host structure
[ht,h1,h2,em,rhom,nim,vol,ef,rhof,nif,alpha]=lay_const(lay,lh,i_em,i_ef,i_rhom,i_rhof,i_nim,i_nif,lvol,lalpha,lmat,lfib);
%bonding layer
if(glueEl)
lhg=coords(nodes(glueEl(1),nx*ny*nz),3)-coords(nodes(glueEl(1),nx*ny),3);
[htg,h1g,h2g,emg,rhomg,nimg,volg,efg,rhofg,nifg,alphag]=lay_const(1,lhg,i_em,i_ef,i_rhom,i_rhof,i_nim,i_nif,0,0,4,1);
end
% titanium alloy rivets - matrix(5)
if(rivetEl)
lhr=coords(nodes(rivetEl(1),nx*ny*nz),3)-coords(nodes(rivetEl(1),nx*ny),3);
[htr,h1r,h2r,emr,rhomr,nimr,volr,efr,rhofr,nifr,alphar]=lay_const(1,lhr,i_em,i_ef,i_rhom,i_rhof,i_nim,i_nif,0,0,5,1);
end
clear i_em i_ef i_rhom i_rhof i_nim i_nif lvol lalpha
% Qpzt=[ 1.2035    0.7518    0.7509         0         0         0
%     0.7518    1.2035    0.7509         0         0         0
%     0.7509    0.7509    1.1087         0         0         0
%          0         0         0    0.2105         0         0
%          0         0         0         0    0.2257         0
%          0         0         0         0         0    0.2105]*1e11;

[ksi,wx]=gll(nx); % weights and nodes distribution
[eta,wy]=gll(ny);
[dzeta,wz]=gll(nz);
% vandermonde approach
[Qx]=Vandermonde_v2(ksi,nx);
[Qy]=Vandermonde_v2(eta,ny);
[Qz]=Vandermonde_v2(dzeta,nz);
%ne=1; % element number
disp('elastic coefficient dependence on evironmental conditions')
D11=zeros(fen*NofElNodes,1);
D12=zeros(fen*NofElNodes,1);
D13=zeros(fen*NofElNodes,1);
D14=zeros(fen*NofElNodes,1);
D22=zeros(fen*NofElNodes,1);
D23=zeros(fen*NofElNodes,1);
D24=zeros(fen*NofElNodes,1);
D33=zeros(fen*NofElNodes,1);
D34=zeros(fen*NofElNodes,1);
D44=zeros(fen*NofElNodes,1);
D55=zeros(fen*NofElNodes,1);
D56=zeros(fen*NofElNodes,1);
D66=zeros(fen*NofElNodes,1);
pzt_bond_num=[pztnum,glueEl,rivetEl];
%hoststrnum=setdiff([1:fen],pztnum);
hoststrnum=setdiff([1:fen],pzt_bond_num);
%hoststrnum=1:fen;
% must be change into parallel mode
I=zeros(fen*NofElNodes,1);

for ne=1:fen
    n1=(ne-1)*NofElNodes+1;
    n2=n1+NofElNodes-1;
    I(n1:n2,1)=nodes(ne,:);
end    
ne=hoststrnum(1);
% host structure
[D,rho]=composite3D_transf(ht,h1,h2,rhom,rhof,em,ef,nim,nif,vol,alpha,lay,Temp(nodes(ne,:)),Moist(nodes(ne,:)),prop,nx,ny,nz,Qz);
D11(:,1)=D(1,1);
D12(:,1)=D(1,2);
D13(:,1)=D(1,3);
D14(:,1)=D(1,4);
D22(:,1)=D(2,2);
D23(:,1)=D(2,3);
D24(:,1)=D(2,4);
D33(:,1)=D(3,3);
D34(:,1)=D(3,4);
D44(:,1)=D(4,4);
D55(:,1)=D(5,5);
D56(:,1)=D(5,6);
D66(:,1)=D(6,6);

% bonding
if(glueEl)
ne=glueEl(1);
[Db,rhob]=composite3D_transf(htg,h1g,h2g,rhomg,rhofg,emg,efg,nimg,nifg,0,0,1,Temp(nodes(ne,:)),Moist(nodes(ne,:)),prop,nx,ny,nz,Qz);
end
% titanium rivet
if(rivetEl)
ne=rivetEl(1);
[Dr,rhor]=composite3D_transf(htr,h1r,h2r,rhomr,rhofr,emr,efr,nimr,nifr,0,0,1,Temp(nodes(ne,:)),Moist(nodes(ne,:)),prop,nx,ny,nz,Qz);
end
clear Temp Moist;
%fen=100;
%% PRELIMINARY CALCULATIONS
disp('preliminary calculations');
mg=zeros(dof,1); % global mass matrix
% global vector of forces induced by the actuator
switch multi
    case 'multipzt'
        [colsensor,cc]=size(pztEl);
        %col= length(outputs);
        col=1;
    case 'multifocus'
        [col,cc]=size(focal);
end
Fa=zeros(dof,col);
%%
disp('assembly mass matrix');
Inc=zeros(fen,3*NofElNodes);
% connectivity matrix Inc by dof
disp('connectivity matrix...');
 for ne=1:fen
     Inc(ne,1:3:end)=3*nodes(ne,:)-2;
     Inc(ne,2:3:end)=3*nodes(ne,:)-1;
     Inc(ne,3:3:end)=3*nodes(ne,:);
 end;

% host structure
for ne=hoststrnum
    %[med]=mass3D_mex(rho,coords(nodes(ne,:),1),coords(nodes(ne,:),2),coords(nodes(ne,:),3),ksi,wx,Qx,nx,dzeta,wz,Qz,nz);
     med=mass3D_v2(rho,coords(nodes(ne,:),1),coords(nodes(ne,:),2),coords(nodes(ne,:),3),Qx,Qy,Qz,ksi,eta,dzeta,wx,wy,wz);
    mg(Inc(ne,:))=mg(Inc(ne,:))+med;% assembly mass matrix
end
% pzt
for ne=pztnum
%    %[mepzt]=mass3D_mex(rho_pzt,coords(nodes(ne,:),1),coords(nodes(ne,:),2),(coords(nodes(ne,:),3)),ksi,wx,Qx,nx,dzeta,wz,Qz,nz);
      mepzt=mass3D_v2(rho_pzt,coords(nodes(ne,:),1),coords(nodes(ne,:),2),coords(nodes(ne,:),3),Qx,Qy,Qz,ksi,eta,dzeta,wx,wy,wz);
    mg(Inc(ne,:))=mg(Inc(ne,:))+mepzt;% assembly mass matrix
end
% bonding layer
for ne=glueEl
%    %[mepzt]=mass3D_mex(rhob,coords(nodes(ne,:),1),coords(nodes(ne,:),2),(coords(nodes(ne,:),3)),ksi,wx,Qx,nx,dzeta,wz,Qz,nz);
      meg=mass3D_v2(rhob,coords(nodes(ne,:),1),coords(nodes(ne,:),2),coords(nodes(ne,:),3),Qx,Qy,Qz,ksi,eta,dzeta,wx,wy,wz);
    mg(Inc(ne,:))=mg(Inc(ne,:))+meg;% assembly mass matrix
end
% rivets
for ne=rivetEl
      meg=mass3D_v2(rhor,coords(nodes(ne,:),1),coords(nodes(ne,:),2),coords(nodes(ne,:),3),Qx,Qy,Qz,ksi,eta,dzeta,wx,wy,wz);
    mg(Inc(ne,:))=mg(Inc(ne,:))+meg;% assembly mass matrix
end
switch mode
    case 'gpu'
    mg=gpuArray(mg); 
    case 'cpu'
        
end
pzt_coords=coords(outputs/3,:);
switch multi
    case 'multifocus'
%% signals for multiple focals
        [ccpzt,cc]=size(pztEl);
        ST=zeros(nft,col,ccpzt);
        dt=tt/nft;   % calculation time step [s]
        for j=1:col
        Dist=sqrt((focal(j,1)-pzt_coords(:,1)).^2+(focal(j,2)-pzt_coords(:,2)).^2+(focal(j,3)-pzt_coords(:,3)).^2);
        % time delay
        t1=(max(Dist)-Dist)/cg;
        %signal    
        for k=1:ccpzt
            [t,st]=Hanning_signal(dt,nft,f_1,f_2,0);
            stt=designed_waveform_fun(st(1:nft/4),Dist(k),dt); % different signal for each pzt actuator
            ST(:,j,k)=stt;
        end
        end
    case 'multipzt'
%% signals for multiple excitations
        ST=zeros(nft,col);
        dt=tt/nft;   % calculation time step [s]
        [t,st]=Hanning_signal(dt,nft,f_1,f_2,0); % the same signal for all pzt actuators
        for k=1:col
        ST(:,k)=st;
        end
end

%  figure;
%  plot(t,st);pause(1);
 
%% forces induced by the pzt actuators

%Fa2=zeros(dof,1);
[Vx,Vzx]=Vandermonde_old(ksi,dzeta,nx,nz);
disp('forces...');
% commented code also works if only one layer of elements is used for pzt
% modelling
% for j=1:length(PZT_actuator)
%     for ne=PZT_actuator(j).pztEl
%         [kufi,kfifi]=coupl3D_v2(epzt,gpzt,coords(nodes(ne,:),1),coords(nodes(ne,:),2),(coords(nodes(ne,:),3)),Qx,Qy,Qz,ksi,eta,dzeta,wx,wy,wz);
%         kfiu=kufi';
%         [fa,fi]=Forces_actuator_mex2(V0,kfiu,kfifi,nx,nz);    
%         Fa(Inc(ne,:))=Fa(Inc(ne,:))+fa; % assembly
%     end
% end
for j=1:length(PZT_actuator)
    ne=PZT_actuator(j).pztEl;
    nodespzt=nodes(ne,:);
%     [p,nodespzt]=remove_free_spec_nodes(coords,nodespzt);
    % electric boundary conditions for sensor in open circuit
    EBC0=nodespzt(:,1:nx*ny);   
    EBCI=nodespzt(:,nx*ny+1:end-nx*ny);   
    EBCN=nodespzt(:,end-nx*ny+1:end);
    EBC0=unique(EBC0);
    EBCI=unique(EBCI);
    EBCN=unique(EBCN);
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
    dofpzt=dof;
    Kfifi=sparse(dofpzt/3,dofpzt/3);
    Kfiu=sparse(dofpzt/3,dofpzt);
    Incpzt=zeros(length(PZT_sensor(j).pztEl),3*NofElNodes);
    
    Incpzt(:,1:3:end)=3*nodespzt(:,:)-2;
    Incpzt(:,2:3:end)=3*nodespzt(:,:)-1;
    Incpzt(:,3:3:end)=3*nodespzt(:,:)-0;
    cpzt=0;
    for ne=PZT_actuator(j).pztEl       
        cpzt=cpzt+1;
        %[kfiu,kfifi]=couplingPZT3D_new_mex2(coords(nodes(ne,:),1),coords(nodes(ne,:),2),(coords(nodes(ne,:),3)),ksi,wx,nx,epzt,gpzt,Vx,dzeta,wz,nz,Vzx);
        [kufi,kfifi]=coupl3D_v2(epzt,gpzt,coords(nodes(ne,:),1),coords(nodes(ne,:),2),(coords(nodes(ne,:),3)),Qx,Qy,Qz,ksi,eta,dzeta,wx,wy,wz);
        kfiu=sparse(kufi');kfifi=sparse(kfifi);
        [iout]=dofs3D(Incpzt,cpzt,dofpzt,NofElNodes);
        [iout2]=dofs3Dfifi(nodespzt,cpzt,dofpzt/3,NofElNodes);
        Kfiu=Kfiu+iout2'*kfiu*iout;
        Kfifi=Kfifi+iout2'*kfifi*iout2;
   end
   KfiuA=[Kfiu(EBCI,BC03),Kfiu(EBCI,BCI3),Kfiu(EBCI,BCN3)];% actuator BC
   KfifiA=Kfifi(EBCI,EBCI); % actuator 
   Kfifiin=Kfifi(EBCI,EBCN);
   V1=sparse(dof/3,1);
   V1(EBCN,1)=V0; % top
   V1(EBC0,1)=0; % bottom
   Fetemp=Kfifiin*V1(EBCN);
   fiAi=-KfifiA\Fetemp; % induced potential
   FiA=sparse(dofpzt/3,1);
   %FiA=[V1(EBC0);fiAi;V1(EBCN)]; % potenitial together with applied voltage
   FiA(EBCN)=V0;
   FiA(EBCI)=fiAi;
   FaL=Kfiu'*FiA; % local forces in pzt actuator j
    PZT_actuator(j).KfiuA=KfiuA;
   PZT_actuator(j).Kfifiin=Kfifiin;
   PZT_actuator(j).KfifiA=KfifiA;
%    PZT_actuator(j).KfifiO=KfifiO;
%    PZT_actuator(j).KfiuO=KfiuO;
%    PZT_actuator(j).KfiuO=KfiuO;
   PZT_actuator(j).Kfiu=Kfiu;
   PZT_actuator(j).EBCI=[EBCI];
   PZT_actuator(j).EBC=[EBCI;EBCN]';
   PZT_actuator(j).EBCN=EBCN;
   PZT_actuator(j).EBC0=EBC0;
   % global numbers
%     ne=PZT_actuator(j).pztEl;
%     EBC0=nodes(ne,1:nx*ny);   
%     EBCI=nodes(ne,nx*ny+1:end-nx*ny);   
%     EBCN=nodes(ne,end-nx*ny+1:end);
%     EBC0=unique(EBC0);
%     EBCI=unique(EBCI);
%     EBCN=unique(EBCN);
%     BC03=zeros(1,3*length(EBC0));
%     BC03(1:3:end)=3*EBC0-2;
%     BC03(2:3:end)=3*EBC0-1;
%     BC03(3:3:end)=3*EBC0-0;
% 
%     BCI3=zeros(1,3*length(EBCI));
%     BCI3(1:3:end)=3*EBCI-2;
%     BCI3(2:3:end)=3*EBCI-1;
%     BCI3(3:3:end)=3*EBCI-0;
% 
%     BCN3=zeros(1,3*length(EBCN));
%     BCN3(1:3:end)=3*EBCN-2;
%     BCN3(2:3:end)=3*EBCN-1;
%     BCN3(3:3:end)=3*EBCN-0;

    BN=[BC03';BCI3';BCN3'];
    PZT_actuator(j).BN=BN;
    Fa(BN,j)=FaL(BN); % global forces from pzt actuator j
end
clear FaL Fetemp
switch mode
    case 'gpu'
    Fa=gpuArray(Fa); 
    case 'cpu'
        
end
%% prepare matrices for sensing
disp('electrical matrices...');
c=0;
for j=1:length(PZT_sensor)
    %[j,length(PZT_sensor)]
    ne=PZT_sensor(j).pztEl;
    nodespzt=nodes(ne,:);
%     [p,nodespzt]=remove_free_spec_nodes(coords,nodespzt);
    % electric boundary conditions for sensor in open circuit
    EBC0=nodespzt(:,1:nx*ny);   
    EBCI=nodespzt(:,nx*ny+1:end-nx*ny);   
    EBCN=nodespzt(:,end-nx*ny+1:end);
    EBC0=unique(EBC0);
    EBCI=unique(EBCI);
    EBCN=unique(EBCN);
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
    dofpzt=3*max(max(nodespzt));
    Kfifi=sparse(dofpzt/3,dofpzt/3);
    Kfiu=sparse(dofpzt/3,dofpzt);
    Incpzt=zeros(length(PZT_sensor(j).pztEl),3*NofElNodes);
    
    Incpzt(:,1:3:end)=3*nodespzt(:,:)-2;
    Incpzt(:,2:3:end)=3*nodespzt(:,:)-1;
    Incpzt(:,3:3:end)=3*nodespzt(:,:)-0;
    cpzt=0;
    for ne=PZT_sensor(j).pztEl       
        cpzt=cpzt+1;
        %[kfiu,kfifi]=couplingPZT3D_new_mex2(coords(nodes(ne,:),1),coords(nodes(ne,:),2),(coords(nodes(ne,:),3)),ksi,wx,nx,epzt,gpzt,Vx,dzeta,wz,nz,Vzx);
        [kufi,kfifi]=coupl3D_v2(epzt,gpzt,coords(nodes(ne,:),1),coords(nodes(ne,:),2),(coords(nodes(ne,:),3)),Qx,Qy,Qz,ksi,eta,dzeta,wx,wy,wz);
        kfiu=sparse(kufi'); kfifi=sparse(kfifi);
        [iout]=dofs3D(Incpzt,cpzt,dofpzt,NofElNodes);
        [iout2]=dofs3Dfifi(nodespzt,cpzt,dofpzt/3,NofElNodes);
        Kfiu=Kfiu+iout2'*kfiu*iout;
        Kfifi=Kfifi+iout2'*kfifi*iout2;
    end
   
   KfifiO=[Kfifi(EBCI,EBCI),Kfifi(EBCI,EBCN);
           Kfifi(EBCN,EBCI),Kfifi(EBCN,EBCN)]; % sensor in open circuit
   KfiuO=[Kfiu(EBCI,BC03),Kfiu(EBCI,BCI3),Kfiu(EBCI,BCN3);
          Kfiu(EBCN,BC03),Kfiu(EBCN,BCI3),Kfiu(EBCN,BCN3)]; 
   
   PZT_sensor(j).Kfifi=Kfifi;
   PZT_sensor(j).Kfiu=Kfiu;
   PZT_sensor(j).KfifiO=KfifiO;
   PZT_sensor(j).KfiuO=KfiuO;
   PZT_sensor(j).Kfiu=Kfiu;
   PZT_sensor(j).EBC=[EBCI;EBCN]';
   PZT_sensor(j).EBCN=EBCN;
   PZT_sensor(j).EBC0=EBC0;
   % global numbers
%    ne=PZT_sensor(j).pztEl;
%     EBC0=nodes(ne,1:nx*ny);   
%     EBCI=nodes(ne,nx*ny+1:end-nx*ny);   
%     EBCN=nodes(ne,end-nx*ny+1:end);
%     EBC0=unique(EBC0);
%     EBCI=unique(EBCI);
%     EBCN=unique(EBCN);
%     BC03=zeros(1,3*length(EBC0));
%     BC03(1:3:end)=3*EBC0-2;
%     BC03(2:3:end)=3*EBC0-1;
%     BC03(3:3:end)=3*EBC0-0;
% 
%     BCI3=zeros(1,3*length(EBCI));
%     BCI3(1:3:end)=3*EBCI-2;
%     BCI3(2:3:end)=3*EBCI-1;
%     BCI3(3:3:end)=3*EBCI-0;
% 
%     BCN3=zeros(1,3*length(EBCN));
%     BCN3(1:3:end)=3*EBCN-2;
%     BCN3(2:3:end)=3*EBCN-1;
%     BCN3(3:3:end)=3*EBCN-0;
     BN=[BC03';BCI3';BCN3'];
    PZT_sensor(j).BN=BN;
end
clear Incpzt Kfifi KfifiA KfifiO Kfifin Kfiu KfiuA KfiuO FiAi FiA;
%% Output file for solution
outfile_voltage=fullfile('outputs',['output',num2str(k_test)],['voltage',num2str(k_test)]);
outfile_displ=fullfile('outputs',['output',num2str(k_test)],['displ',num2str(k_test)]);
outfile_time=fullfile('outputs',['output',num2str(k_test)],['time',num2str(k_test)]);
% format for output
% str1=[];
% for kk=1:length(PZT_sensor);
%     str1=[str1,'%17.9e '];
% end
% str1=[str1,'\n'];

%%
disp('local derivatives');
[Nprimx,Nprimy,Nprimz]=shape3D_prim_v2p(nx,ny,nz,Qx,Qy,Qz,ksi',eta',dzeta');


[indxi,indxj]=find(Nprimx);% indices of nonzero values for derivative in x direction
[indyi,indyj]=find(Nprimy);% indices of nonzero values for derivative in y direction
[indzi,indzj]=find(Nprimz);% indices of nonzero values for derivative in z direction
ixl=length(indxi);
iyl=length(indyi);
izl=length(indzi);
Npxs=zeros(ixl,1);
Npys=zeros(iyl,1);
Npzs=zeros(izl,1);
for k=1:ixl
    Npxs(k)=Nprimx(indxi(k),indxj(k));      
end
for k=1:iyl
    Npys(k)=Nprimy(indyi(k),indyj(k));
end
for k=1:izl
    Npzs(k)=Nprimz(indzi(k),indzj(k));
end
Nxs=repmat(Npxs,fen,1);
Nys=repmat(Npys,fen,1);
Nzs=repmat(Npzs,fen,1);
clear Npxs Npys Npzs 
Indxi=zeros(fen*ixl,1);
Indxj=zeros(fen*ixl,1);
Indyi=zeros(fen*iyl,1);
Indyj=zeros(fen*iyl,1);
Indzi=zeros(fen*izl,1);
Indzj=zeros(fen*izl,1);
for ne=1:fen
    offset=(ne-1)*NofElNodes;
    n1=(ne-1)*ixl+1;
    n2=n1+ixl-1;
    n3=(ne-1)*iyl+1;
    n4=n3+iyl-1;
    n5=(ne-1)*izl+1;
    n6=n5+izl-1;
    Indxi(n1:n2,1)=offset+indxi;
    Indxj(n1:n2,1)=offset+indxj;
    Indyi(n3:n4,1)=offset+indyi;
    Indyj(n3:n4,1)=offset+indyj;
    Indzi(n5:n6,1)=offset+indzi;
    Indzj(n5:n6,1)=offset+indzj;
end
Npx=sparse(Indxi,Indxj,Nxs);
clear Nxs;
Npy=sparse(Indyi,Indyj,Nys);
clear Nys;
Npz=sparse(Indzi,Indzj,Nzs);
clear Nzs;
clear indxi indxj indyi indyj indzi indzj;
X=zeros(fen*NofElNodes,1);
Y=zeros(fen*NofElNodes,1);
Z=zeros(fen*NofElNodes,1);
for ne=1:fen
    n1=(ne-1)*NofElNodes+1;
    n2=n1+NofElNodes-1;
    X(n1:n2)=coords(nodes(ne,:),1);
    Y(n1:n2)=coords(nodes(ne,:),2);
    Z(n1:n2)=coords(nodes(ne,:),3);
end
switch mode
    case 'gpu'
    Npx=gpuArray(Npx);
    Npy=gpuArray(Npy);
    Npz=gpuArray(Npz);
    
    X=gpuArray(X);
    Y=gpuArray(Y);
    Z=gpuArray(Z);
   
    case 'cpu'
        
end
disp('jacobians'); 
J11=Npx*X;
J21=Npx*Y;
J31=Npx*Z;
J12=Npy*X;
J22=Npy*Y;
J32=Npy*Z;
J13=Npz*X;
J23=Npz*Y;
J33=Npz*Z;

y=[coords(nodes(1,2),1),coords(nodes(1,2),2),coords(nodes(1,2),3)];
x=[coords(nodes(1,1),1),coords(nodes(1,1),2),coords(nodes(1,1),3)];
h=sum((x-y).^2).^0.5;
clear X Y Z
clear Indxi Indxj Indyi Indyj Indzi Indzj
%disp('...done');    


disp('elastic constants');


% overwrite mechanical properties in pzt elements
for ne=pztnum
    n1=(ne-1)*NofElNodes+1;
     n2=n1+NofElNodes-1;
     D11(n1:n2)=Qpzt(1,1);
     D12(n1:n2)=Qpzt(1,2);
     D13(n1:n2)=Qpzt(1,3);
     D14(n1:n2)=Qpzt(1,4);
     D22(n1:n2)=Qpzt(2,2);
     D23(n1:n2)=Qpzt(2,3);
     D24(n1:n2)=Qpzt(2,4);
     D33(n1:n2)=Qpzt(3,3);
     D34(n1:n2)=Qpzt(3,4);
     D44(n1:n2)=Qpzt(4,4);
     D55(n1:n2)=Qpzt(5,5);
     D56(n1:n2)=Qpzt(5,6);
     D66(n1:n2)=Qpzt(6,6);
end
% bonding layer
% overwrite mechanical properties in bonding layer elements
for ne=glueEl
    n1=(ne-1)*NofElNodes+1;
    n2=n1+NofElNodes-1;
    D11(n1:n2,1)=Db(1,1);
    D12(n1:n2,1)=Db(1,2);
    D13(n1:n2,1)=Db(1,3);
    D14(n1:n2,1)=Db(1,4);
    D22(n1:n2,1)=Db(2,2);
    D23(n1:n2,1)=Db(2,3);
    D24(n1:n2,1)=Db(2,4);
    D33(n1:n2,1)=Db(3,3);
    D34(n1:n2,1)=Db(3,4);
    D44(n1:n2,1)=Db(4,4);
    D55(n1:n2,1)=Db(5,5);
    D56(n1:n2,1)=Db(5,6);
    D66(n1:n2,1)=Db(6,6);
end
% rivets
% overwrite mechanical properties in rivet elements
for ne=rivetEl
    n1=(ne-1)*NofElNodes+1;
    n2=n1+NofElNodes-1;
    D11(n1:n2,1)=Dr(1,1);
    D12(n1:n2,1)=Dr(1,2);
    D13(n1:n2,1)=Dr(1,3);
    D14(n1:n2,1)=Dr(1,4);
    D22(n1:n2,1)=Dr(2,2);
    D23(n1:n2,1)=Dr(2,3);
    D24(n1:n2,1)=Dr(2,4);
    D33(n1:n2,1)=Dr(3,3);
    D34(n1:n2,1)=Dr(3,4);
    D44(n1:n2,1)=Dr(4,4);
    D55(n1:n2,1)=Dr(5,5);
    D56(n1:n2,1)=Dr(5,6);
    D66(n1:n2,1)=Dr(6,6);
end

%disp('done');
disp('determinant and inverse Jacobian')
%%%
% /* determinant and inverse Jacobian  */
[invJ11,invJ12,invJ13,invJ21,invJ22,invJ23,invJ31,invJ32,invJ33]=inv_jacp(J11,J12,J13,J21,J22,J23,J31,J32,J33);
detJ=det_jacp(J11,J12,J13,J21,J22,J23,J31,J32,J33);
switch mode
    case 'gpu'
        j11=gather(J11);j21=gather(J21);j12=gather(J12);j22=gather(J22);j13=gather(J13);j23=gather(J23);j31=gather(J31);j32=gather(J31);j33=gather(J33);
        invj11=gather(invJ11);invj21=gather(invJ21);invj12=gather(invJ12);invj22=gather(invJ22);
        invj13=gather(invJ13);invj23=gather(invJ23);invj31=gather(invJ31);invj32=gather(invJ31);invj33=gather(invJ33);
        det_jac=detJ;
        save([meshfile(1:end-4),'_jacobians'],'j11','j21','j12','j22','j13','j23','j31','j32','j33','invj11','invj12','invj21','invj22','invj13','invj23','invj31','invj32','invj33','det_jac');
    case 'cpu'
        j11=J11;j21=J21;j12=J12;j22=J22;j13=J13;j23=J23;j31=J31;j32=J32;j33=J33;
        invj11=invJ11;invj21=invJ21;invj12=invJ12;invj22=invJ22;
        invj13=invJ13;invj23=invJ23;invj31=invJ31;invj32=invJ32;invj33=invJ33;
        save([meshfile(1:end-4),'_jacobians'],'j11','j21','j12','j22','j13','j23','j31','j32','j33','invj11','invj12','invj21','invj22','invj13','invj23','invj31','invj32','invj33','det_jac');
end
clear j11 j12 j21 j22 j13 j23 j32 j33 invj11 invj12 invj21 invj22 invj13 invj23 invj32 invj33 det_jac
%disp('done');
%
D11=repmat(D11,1,col);
D12=repmat(D12,1,col);
D13=repmat(D13,1,col);
D14=repmat(D14,1,col);
D22=repmat(D22,1,col);
D23=repmat(D23,1,col);
D24=repmat(D24,1,col);
D33=repmat(D33,1,col);
D34=repmat(D34,1,col);
D44=repmat(D44,1,col);
D55=repmat(D55,1,col);
D56=repmat(D56,1,col);
D66=repmat(D66,1,col);
invJ11=repmat(invJ11,1,col);
invJ12=repmat(invJ12,1,col);
invJ13=repmat(invJ13,1,col);
invJ21=repmat(invJ21,1,col);
invJ22=repmat(invJ22,1,col);
invJ23=repmat(invJ23,1,col);
invJ31=repmat(invJ31,1,col);
invJ32=repmat(invJ32,1,col);
invJ33=repmat(invJ33,1,col);
switch mode
    case 'gpu'
    D11=gpuArray(D11);D12=gpuArray(D12);D13=gpuArray(D13);D14=gpuArray(D14);
    D22=gpuArray(D22);D23=gpuArray(D23);D24=gpuArray(D24);
    D33=gpuArray(D33);D34=gpuArray(D34);
    D44=gpuArray(D44);D55=gpuArray(D55);D56=gpuArray(D56);D66=gpuArray(D66);
    
    case 'cpu'
        
end
%%


%%
      
% weights at quadrature points
cc=0;
www=zeros(NofElNodes,1);
for kz=1:nz % dzeta
    for ky=1:ny % eta
        for kx=1:nx % ksi  
            cc=cc+1;
            www(cc,1)=wx(kx)*wy(ky)*wz(kz);
        end
    end
end
switch mode
    case 'gpu'
    www=gpuArray(www);
    
    case 'cpu'
        
end
% 
WWW=repmat(www,fen,1);
WWWDetJ=WWW.*detJ;
WWWDetJ=repmat(WWWDetJ,1,col);

clear WWW www detJ
clear J11 J12 J13 J21 J22 J23 J31 J32 J33
% global 
% WWWDetJglob=zeros(dof,1);
% for ne=1:fen
%     n1=(ne-1)*NofElNodes+1;
%     n2=n1+NofElNodes-1;
%     WWWDetJglob(nodes(ne,:))=WWWDetJ(n1:n2);
% end
% end of global
disp('preliminary calculations');       
Iu=zeros(fen*NofElNodes,1);
Iv=zeros(fen*NofElNodes,1);
Iw=zeros(fen*NofElNodes,1);
for ne=1:fen
    n1=(ne-1)*NofElNodes+1;
    n2=n1+NofElNodes-1;
    Iu(n1:n2,1)=3*nodes(ne,:)-2;
    Iv(n1:n2,1)=3*nodes(ne,:)-1;
    Iw(n1:n2,1)=3*nodes(ne,:);
end    

Uc=zeros(dof,col); % displacement on cpu
switch mode
    case 'gpu'
    Iu=gpuArray(uint32(Iu));
    Iv=gpuArray(uint32(Iv));
    Iw=gpuArray(uint32(Iw));       
    U=zeros(dof,col,'double','gpuArray');
    UX=zeros(fen*NofElNodes,col,'double','gpuArray');
    UY=zeros(fen*NofElNodes,col,'double','gpuArray');
    UZ=zeros(fen*NofElNodes,col,'double','gpuArray');
    FF=zeros(3*fen*NofElNodes,col,'gpuArray');
    
    case 'cpu'
    U=zeros(dof,col);
    UX=zeros(fen*NofElNodes,col);
    UY=zeros(fen*NofElNodes,col);
    UZ=zeros(fen*NofElNodes,col);  
    FF=zeros(3*fen*NofElNodes,col);
end

%% Pre Main       
[I_G_dofs,I_L_dofs] = global_local_12baskets_3dofs_v2(I_G,I_L,mode);

% Itemp=zeros(3*length(IG1),1);
% Itemp(1:3:end)=3*IG1-2;
% Itemp(2:3:end)=3*IG1-1;
% Itemp(3:3:end)=3*IG1;
% IG1=Itemp;
% 
% Itemp=zeros(3*length(IG2),1);
% Itemp(1:3:end)=3*IG2-2;
% Itemp(2:3:end)=3*IG2-1;
% Itemp(3:3:end)=3*IG2;
% IG2=Itemp;
% 
% Itemp=zeros(3*length(IG3),1);
% Itemp(1:3:end)=3*IG3-2;
% Itemp(2:3:end)=3*IG3-1;
% Itemp(3:3:end)=3*IG3;
% IG3=Itemp;
% 
% Itemp=zeros(3*length(IG4),1);
% Itemp(1:3:end)=3*IG4-2;
% Itemp(2:3:end)=3*IG4-1;
% Itemp(3:3:end)=3*IG4;
% IG4=Itemp;
% 
% Itemp=zeros(3*length(IG5),1);
% Itemp(1:3:end)=3*IG5-2;
% Itemp(2:3:end)=3*IG5-1;
% Itemp(3:3:end)=3*IG5;
% IG5=Itemp;
% 
% Itemp=zeros(3*length(IG6),1);
% Itemp(1:3:end)=3*IG6-2;
% Itemp(2:3:end)=3*IG6-1;
% Itemp(3:3:end)=3*IG6;
% IG6=Itemp;
% 
% Itemp=zeros(3*length(IG7),1);
% Itemp(1:3:end)=3*IG7-2;
% Itemp(2:3:end)=3*IG7-1;
% Itemp(3:3:end)=3*IG7;
% IG7=Itemp;
% 
% Itemp=zeros(3*length(IG8),1);
% Itemp(1:3:end)=3*IG8-2;
% Itemp(2:3:end)=3*IG8-1;
% Itemp(3:3:end)=3*IG8;
% IG8=Itemp;
% 
% Itemp=zeros(3*length(IG9),1);
% Itemp(1:3:end)=3*IG9-2;
% Itemp(2:3:end)=3*IG9-1;
% Itemp(3:3:end)=3*IG9;
% IG9=Itemp;
% 
% Itemp=zeros(3*length(IG10),1);
% Itemp(1:3:end)=3*IG10-2;
% Itemp(2:3:end)=3*IG10-1;
% Itemp(3:3:end)=3*IG10;
% IG10=Itemp;
% 
% Itemp=zeros(3*length(IG11),1);
% Itemp(1:3:end)=3*IG11-2;
% Itemp(2:3:end)=3*IG11-1;
% Itemp(3:3:end)=3*IG11;
% IG11=Itemp;
% 
% Itemp=zeros(3*length(IG12),1);
% Itemp(1:3:end)=3*IG12-2;
% Itemp(2:3:end)=3*IG12-1;
% Itemp(3:3:end)=3*IG12;
% IG12=Itemp;
% 
% Itemp=zeros(3*length(IL1),1);
% Itemp(1:3:end)=3*IL1-2;
% Itemp(2:3:end)=3*IL1-1;
% Itemp(3:3:end)=3*IL1;
% IL1=Itemp;
% 
% Itemp=zeros(3*length(IL2),1);
% Itemp(1:3:end)=3*IL2-2;
% Itemp(2:3:end)=3*IL2-1;
% Itemp(3:3:end)=3*IL2;
% IL2=Itemp;
% 
% 
% Itemp=zeros(3*length(IL3),1);
% Itemp(1:3:end)=3*IL3-2;
% Itemp(2:3:end)=3*IL3-1;
% Itemp(3:3:end)=3*IL3;
% IL3=Itemp;
% 
% Itemp=zeros(3*length(IL4),1);
% Itemp(1:3:end)=3*IL4-2;
% Itemp(2:3:end)=3*IL4-1;
% Itemp(3:3:end)=3*IL4;
% IL4=Itemp;
% 
% Itemp=zeros(3*length(IL5),1);
% Itemp(1:3:end)=3*IL5-2;
% Itemp(2:3:end)=3*IL5-1;
% Itemp(3:3:end)=3*IL5;
% IL5=Itemp;
% 
% Itemp=zeros(3*length(IL6),1);
% Itemp(1:3:end)=3*IL6-2;
% Itemp(2:3:end)=3*IL6-1;
% Itemp(3:3:end)=3*IL6;
% IL6=Itemp;
% 
% Itemp=zeros(3*length(IL7),1);
% Itemp(1:3:end)=3*IL7-2;
% Itemp(2:3:end)=3*IL7-1;
% Itemp(3:3:end)=3*IL7;
% IL7=Itemp;
% 
% Itemp=zeros(3*length(IL8),1);
% Itemp(1:3:end)=3*IL8-2;
% Itemp(2:3:end)=3*IL8-1;
% Itemp(3:3:end)=3*IL8;
% IL8=Itemp;
% 
% Itemp=zeros(3*length(IL9),1);
% Itemp(1:3:end)=3*IL9-2;
% Itemp(2:3:end)=3*IL9-1;
% Itemp(3:3:end)=3*IL9;
% IL9=Itemp;
% 
% Itemp=zeros(3*length(IL10),1);
% Itemp(1:3:end)=3*IL10-2;
% Itemp(2:3:end)=3*IL10-1;
% Itemp(3:3:end)=3*IL10;
% IL10=Itemp;
% 
% Itemp=zeros(3*length(IL11),1);
% Itemp(1:3:end)=3*IL11-2;
% Itemp(2:3:end)=3*IL11-1;
% Itemp(3:3:end)=3*IL11;
% IL11=Itemp;
% 
% Itemp=zeros(3*length(IL12),1);
% Itemp(1:3:end)=3*IL12-2;
% Itemp(2:3:end)=3*IL12-1;
% Itemp(3:3:end)=3*IL12;
% IL12=Itemp;
% 
% clear Itemp;
% 
% switch mode
%     case 'gpu'
%     IG1=gpuArray(IG1);
%     IG2=gpuArray(IG2);
%     IG3=gpuArray(IG3);
%     IG4=gpuArray(IG4);
%     IG5=gpuArray(IG5);
%     IG6=gpuArray(IG6);
%     IG7=gpuArray(IG7);
%     IG8=gpuArray(IG8);
%     IG9=gpuArray(IG9);
%     IG10=gpuArray(IG10);
%     IG11=gpuArray(IG11);
%     IG12=gpuArray(IG12);
% 
%     IL1=gpuArray(IL1);
%     IL2=gpuArray(IL2);
%     IL3=gpuArray(IL3);
%     IL4=gpuArray(IL4);
%     IL5=gpuArray(IL5);
%     IL6=gpuArray(IL6);
%     IL7=gpuArray(IL7);
%     IL8=gpuArray(IL8);
%     IL9=gpuArray(IL9);
%     IL10=gpuArray(IL10);
%     IL11=gpuArray(IL11);
%     IL12=gpuArray(IL12);
% 
%     case 'cpu'
%     
% end

%% Boundary conditions ??
%  apply to Npx, Npy, Npz ?? Npxt
%%
Npxt=Npx';
Npyt=Npy';
Npzt=Npz';

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

dv(1:2:end)=etad_xy; % in-plane damping factor
dv(2:2:end)=etad_xy; % in-plane damping factor
dv(3:3:end)=etad_z; % out-of-plane damping factor
cg = dv.*mg; % damping matrix
clear dv;
switch mode
    case 'gpu'
    uold=zeros(dof,col,'double','gpuArray');
    unew=zeros(dof,col,'double','gpuArray'); % memory prealocation
    v=zeros(dof,col,'double','gpuArray');
    case 'cpu'
    uold=zeros(dof,col);
    unew=zeros(dof,col); % memory prealocation
    v=zeros(dof,col);
end


mg2=a2.*mg;
mg=a0.*mg;
mg0=1./(mg+a1*cg);

mg2=repmat(mg2,1,col);
mg=repmat(mg,1,col);
mg0=repmat(mg0,1,col);
cg=repmat(cg,1,col);

voltage=zeros(nft,length(PZT_sensor)*length(PZT_actuator));
displ=zeros(nft,length(outputs));
%clear coords D V1

switch mode
    case 'gpu'
    nodes=gpuArray(uint32(nodes));
    ne=zeros(1,1,'uint32','gpuArray'); 
    NofElNodes=gpuArray(NofElNodes);
    n1=zeros(1,1,'uint32','gpuArray');
    n2=zeros(1,1,'uint32','gpuArray');
    case 'cpu'
%     IG=[IG1;IG2;IG3;IG4;IG5;IG6;IG7;IG8;IG9;IG10;IG11;IG12];
%     IL=[IL1;IL2;IL3;IL4;IL5;IL6;IL7;IL8;IL9;IL10;IL11;IL12];
end
%
 switch mode
        case 'gpu'
            F1=zeros(dof,col,'double','gpuArray');  
        case 'cpu'
            F1=zeros(dof,col);
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

% transform from global displacement matrix to local elemental metrix (multiple right hand sides)
    for k=1:col
        UX(:,k)=U(Iu,k);  
        UY(:,k)=U(Iv,k);
        UZ(:,k)=U(Iw,k);
    end
    
   
    % calculate strains
    [Exx,Eyy,Ezz,Exy,Eyz,Exz]=strains_spec_p(Npx,Npy,Npz,invJ11,invJ12,invJ13,invJ21,invJ22,invJ23,invJ31,invJ32,invJ33,UX,UY,UZ);

    
    
    % calculatate stresses                 
    [Sxx,Syy,Szz,Sxy,Syz,Sxz]=stresses_spec_p(Exx,Eyy,Ezz,Exy,Eyz,Exz,D11,D12,D13,D14,D22,D23,D24,D33,D34,D44,D55,D56,D66);
   
    clear Exx Eyy Ezz Exy Eyz Exz;
    % calculatate forces
    fu=Npxt*(Sxx.*invJ11.*WWWDetJ)+Npyt*(Sxx.*invJ21.*WWWDetJ)+Npzt*(Sxx.*invJ31.*WWWDetJ);
    fu=fu+Npxt*(Sxy.*invJ12.*WWWDetJ)+Npyt*(Sxy.*invJ22.*WWWDetJ)+Npzt*(Sxy.*invJ32.*WWWDetJ);
    fu=fu+Npxt*(Sxz.*invJ13.*WWWDetJ)+Npyt*(Sxz.*invJ23.*WWWDetJ)+Npzt*(Sxz.*invJ33.*WWWDetJ);
    
    fv=Npxt*(Syy.*invJ12.*WWWDetJ)+Npyt*(Syy.*invJ22.*WWWDetJ)+Npzt*(Syy.*invJ32.*WWWDetJ);
    fv=fv+Npxt*(Sxy.*invJ11.*WWWDetJ)+Npyt*(Sxy.*invJ21.*WWWDetJ)+Npzt*(Sxy.*invJ31.*WWWDetJ);
    fv=fv+Npxt*(Syz.*invJ13.*WWWDetJ)+Npyt*(Syz.*invJ23.*WWWDetJ)+Npzt*(Syz.*invJ33.*WWWDetJ);
    
    fw=Npxt*(Szz.*invJ13.*WWWDetJ)+Npyt*(Szz.*invJ23.*WWWDetJ)+Npzt*(Szz.*invJ33.*WWWDetJ);
    fw=fw+Npxt*(Syz.*invJ12.*WWWDetJ)+Npyt*(Syz.*invJ22.*WWWDetJ)+Npzt*(Syz.*invJ32.*WWWDetJ);
    fw=fw+Npxt*(Sxz.*invJ11.*WWWDetJ)+Npyt*(Sxz.*invJ21.*WWWDetJ)+Npzt*(Sxz.*invJ31.*WWWDetJ);
    
    clear Sxx Syy Szz Sxy Syz Sxz;
 
    switch mode
        case 'gpu'
        F=zeros(dof,col,'double','gpuArray');
        
        case 'cpu'
        F=zeros(dof,col);
    end
    
    
 switch assembly
     case 'trad'
          
            for ne=1:fen
                  n1=(ne-1)*NofElNodes+1;
                  n2=n1+NofElNodes-1;              
                  F(3*nodes(ne,:)-2,:)=F(3*nodes(ne,:)-2,:)+fu(n1:n2,:);
                  F(3*nodes(ne,:)-1,:)=F(3*nodes(ne,:)-1,:)+fv(n1:n2,:);
                  F(3*nodes(ne,:),:)=F(3*nodes(ne,:),:)+fw(n1:n2,:);
             end 
        clear fu fv fw
       
     case 'mesh'
         
    FF(1:3:end,:)=fu;
    FF(2:3:end,:)=fv;
    FF(3:3:end,:)=fw;
    clear fu fv fw
    for ibasket=1:12
         F(I_G_dofs(:,ibasket))=F(I_G_dofs(:,ibasket))+FF(I_L_dofs(:,ibasket)); 
    end
    %
%     F(IG1,:)=F(IG1,:)+FF(IL1,:);
%     
%     %
%     F(IG2,:)=F(IG2,:)+FF(IL2,:);
%     
%     F(IG3,:)=F(IG3,:)+FF(IL3,:);
%     
%     F(IG4,:)=F(IG4,:)+FF(IL4,:);
%     
%     F(IG5,:)=F(IG5,:)+FF(IL5,:);
%     
%     F(IG6,:)=F(IG6,:)+FF(IL6,:);
%     
%     F(IG7,:)=F(IG7,:)+FF(IL7,:);
%     
%     F(IG8,:)=F(IG8,:)+FF(IL8,:);
%     
%     F(IG9,:)=F(IG9,:)+FF(IL9,:);
%     
%     F(IG10,:)=F(IG10,:)+FF(IL10,:);
%     
%     F(IG11,:)=F(IG11,:)+FF(IL11,:);
%     
%     F(IG12,:)=F(IG12,:)+FF(IL12,:);
    
 end
    %% time integration
    %F3=-F+mg2.*U+Fa*st(nn)-mg.*uold-cg.*v;
    switch multi
        case 'multipzt'
            for j=1:col
                %
                %F1(outputs(j),j)=ST(nn,j)*V0; % force excitation
                pztdof=PZT_actuator(j).BN;
                F1(pztdof,j)=Fa(pztdof,j)*ST(nn,j);
            end
        case 'multifocus'
            for k=1:col
                for j=1:ccpzt
                    pztdof=PZT_actuator(j).BN;
                    F1(pztdof,k)=Fa(pztdof,k)*ST(nn,k,j);
                end
            end
    end
    %F=-F+mg2.*U+F1-mg.*uold;
    F=-F+mg2.*U+F1-mg.*uold+a1*cg.*uold; % piezoelectric excitation with damping

    %F=-F+mg2.*U+F1*st1(nn)+F2*st2(nn)+F3*st3(nn)+F4*st4(nn)-mg.*uold;
    
    %
    unew=mg0.*F;
    
    %
    %%%%
    v = a1*(-uold+unew);  % update velocity
    uold=U;  
    
    %;
    U = unew; 
    
     %% output signals
     switch multi
        case 'multipzt'
            id=0;
            for i1=1:col
                for j1=1:length(outputs)
                id=id+1;
                    switch mode
                        case 'gpu'
                        uc=gather(U(PZT_sensor(j1).BN,i1));
                        uc1=gather(U(outputs(j1),i1));
                        case 'cpu'
                        uc=U(PZT_sensor(j1).BN,i1);
                        uc1=U(outputs(j1),i1);
                    end
                    Ftmp=PZT_sensor(j1).KfiuO*uc;
                    fio=-PZT_sensor(j1).KfifiO\Ftmp;
                    FiS=sparse(dof/3,1);
                    FiS( PZT_sensor(j1).EBC)=fio;
                    volt=mean(fio(end-nx*ny+1:end));
                    voltage(nn,id)=volt;
                    displ(nn,id)=uc1;
                end
            end
        case 'multifocus'
            for j=1:col
                switch mode
                        case 'gpu'
                        uc1=gather(U(outputs,j));
                        case 'cpu'
                        uc1=U(outputs,j);
                end
                displ(nn,(j-1)*ccpzt+1:j*ccpzt)=uc1;
            end
     end
    %% output frame
    % save frame to file
    if (mod(nn,frm_int) == 0)
        c=c+1;
        t_frames(c)=t(nn);
        Uc=gather(U);
         if( max(Uc)> 10*h(1)) disp('integration error'); 
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
        
        outfile_Ux=fullfile('outputs',['output',num2str(k_test)],['Ux_frame' ,num2str(nn)]);
        outfile_Uy=fullfile('outputs',['output',num2str(k_test)],['Uy_frame' ,num2str(nn)]);
        outfile_Uz=fullfile('outputs',['output',num2str(k_test)],['Uz_frame' ,num2str(nn)]);
        
        %save(outfile,'Uc','-ascii'); % frame output  
        Ux = Uc(1:3:end); Uy = Uc(2:3:end);Uz = Uc(3:3:end);
        save(outfile_Ux,'Ux'); % frame output for global vector of Ux displacement
        save(outfile_Uy,'Uy'); % frame output for global vector of Uy displacement
        save(outfile_Uz,'Uz'); % frame output for global vector of Uz displacement
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
save(outfile_voltage,'voltage');   
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




