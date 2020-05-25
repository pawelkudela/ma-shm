disp('.. Reading input data');
%% Signal definition
nft=2^16;    % total number of samples
tt=0.4e-3;   % total calculation time [s] % 
t_1=0e-4;    % excitation initiation time [s]
f_1=100e3/5;   % frequency of the modulation signal [Hz]
f_2=5*f_1;   % frequency of the carrier signal [Hz]
frames=512; % number of frames for movie
frm_int=floor(nft/(frames)); % save displacement with interval time step frm_int (frames)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% properties of composite material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix(1) - aluminium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Young modulus [Pa]';                      i_em(1) = 71e9;
'Poisson ratio';                           i_nim(1) = 0.33;
'density [kg/m3]';                         i_rhom(1) = 2700;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix(2) - A36 steel http://en.wikipedia.org/wiki/A36_steel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Young modulus [Pa]';                      i_em(2) = 200e9;
'Poisson ratio';                           i_nim(2) = 0.26;
'density [kg/m3]';                         i_rhom(2) = 7800;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix(3) - epoxy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Young modulus [Pa]';                      i_em(3) = 3.43e9;
'Poisson ratio';                           i_nim(3) = 0.35;
'density [kg/m3]';                         i_rhom(3) = 1250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix(4) - pzt bonding layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Young modulus [Pa]';                      i_em(4) = 2e9;
'Poisson ratio';                           i_nim(4) = 0.34;
'density [kg/m3]';                         i_rhom(4) = 1170;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix(5) - titanium alloy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Young modulus [Pa]';                      i_em(5) = 110e9;
'Poisson ratio';                           i_nim(5) = 0.34;
'density [kg/m3]';                         i_rhom(5) = 4430;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fibres(1) - fiber glass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Young modulus [Pa]';                      i_ef(1) = 85e9; 
'Poisson ratio';                           i_nif(1) = 0.23; 
'density [kg/m3]';                         i_rhof(1) = 2250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fibres(2) - fiber graphite/carbon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Young modulus [Pa]';                      i_ef(2) = 275.6e9; 
'Poisson ratio';                           i_nif(2) = 0.2; 
'density [kg/m3]';                         i_rhof(2) = 1900;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% properites of material layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'number of material layers';               lay = 4;
'stack thickness sequence [mm]';           lh = [1.5 1.5 1.5 1.5]/lay/1000;
'stack angle sequence [deg]';              lalpha = [0 0 0 0];
'stack matrix sequnce';                    lmat = [3 3 3 3];
'stack fibres sequece';                    lfib = [2 2 2 2];
'volume fraction of fibres';               lvol = [0.5 0.5 0.5 0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% damping coefficient proportional to mass matrix
etad_xy=2e4; % damping ratio
etad_z=2e4; % damping ratio
% geometry definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V0=70; %[V] peak voltage 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% properties of piezoelectric material
% compliance of the pzt when the electric field is constant
% convention Ex,Ey,Ez,Gxz,Gyz,Gxy
% http://www.americanpiezo.com/apc-materials/piezoelectric-properties.html
% Navy I 840
 YE11=8e10; %[N/m^2] youngs modulus at constant electric field
 YE33=6.8e10; %[N/m^2] youngs modulus at constant electric field
 nu11=0.31; % Poisson ratio - assumed

 Spzt=[ 1/YE11      -nu11/YE11    -nu11/YE33   0           0        0;
       -nu11/YE11   1/YE11      -nu11/YE33     0           0        0;
       -nu11/YE33    -nu11/YE33    1/YE33      0           0        0;
        0              0           0       (1+nu11)/YE33   0        0;
        0              0           0           0     (1+nu11)/YE33   0;
        0              0           0           0         0    (1+nu11)/YE11];%[m^2/N] efunda

% matrix of piezoelectric coupling constants - charge constants
% http://www.americanpiezo.com/apc-materials/piezoelectric-properties.html
% Navy I 840
      %sx    sy  sz  sxz syz sxy
dp=  [ 0     0   0   0   480 0;
       0     0   0   480 0   0;
      -175 -175  290 0   0   0]*10^-12;%[C/N] or [m/V]
% voltage constants
      %sx    sy  sz  sxz syz sxy
gp=  [ 0     0   0   0   38 0;
       0     0   0   38 0   0;
      -12.4 -12.4  26.5 0   0   0]*10^-3;%[m^2/C] or [Vm/N] 
% gpzt must be calculated   ???
% permittivity matrix (strain-charge form)
epsT=   [dp(1,5)/gp(1,5)    0                  0;
        0              dp(1,5)/gp(1,5)       0;
        0                 0              dp(3,3)/gp(3,3)]; %[F/m]
% density
rho_pzt=7600;%[kg/m3]
%% boundary conditions
% global b.c. vector with blocked degrees of freedom
BC=[];
%% mesh
meshfile='mesh\sensor_opt_pzt_stiffener_low_rivets_delam_3D.mat'; % 
coords = coords3D;
nodes = nodes3D;
% meshfile contains two matrices: coords and nodes and list of vectors with
% indices in local and global level
% element type: number of nodes in x, y and z direction respectively
% must be compatible with the mesh provided !
nx=6;ny=6;nz=3;
%% environmental effects
envfile='mesh/environment.mat';
%% through thickness integration method
% options: smeared or layered
prop='smeared'; 
%prop='layered';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actuators and sensors: structure defines which element numbers belongs to
% specific actuator and sensor
cd ..;
load('mesh\pztnum_sensor_opt_pzt_stiffener_low_1lay_no_glue_no_stiffener_3D.mat'); %
pztnum=reshape(pztEl,1,[]);
%pztnum=[];
c=0;
PZT_actuator=[];
PZT_sensor=[];
for ne=3:3
   c=c+1;
   PZT_actuator(c).pztEl=[pztEl(ne,:)]; % list of element numbers for pzt no c
end
c=0;
for ne=1:8
   c=c+1;
   PZT_sensor(c).pztEl=[pztEl(c,:)]; % 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bonding layer
load('mesh\gluenum_sensor_opt_pzt_stiffener_low_1lay_no_glue_no_stiffener_3D.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outputs
load('mesh\outputs_sensor_opt_pzt_stiffener_low_1lay_no_glue_no_stiffener_3D.mat');
% % stiffener
% load('mesh\stiffenernum_sensor_opt_pzt_stiffener_low_1lay_no_glue_no_stiffener_3D.mat');