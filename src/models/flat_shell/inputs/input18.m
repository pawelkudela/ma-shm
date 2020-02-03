disp('.. Reading input data');
%% Signal definition
%nft=1024*64;      % total number of samples
nft=150000;      % total number of samples
tt=0.75e-3;          % total calculation time [s] % 
t_1=0e-4;           % excitation initiation time [s]
f_1=50e3/5;        % frequency of the modulation signal [Hz]
f_2=5*f_1;          % frequency of the carrier signal [Hz]
nFrames=512;        % number of frames for animation
field_variable = 'velocity'; % field_variable for saving output data, string: 'displacement', 'velocity', 'acceleration' or 'all'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% delamination
% can be defined inside meshfile
delamination_layer = [3]; % delamination layer counting from the top after which the nodes are splitted (another column for second delamination)
damage_degree = 0.5;  % damage_degree - double in range 0:1 (0 no damage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% properties of composite material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C=[52.55 6.51    5.94     0       0     0;
    6.51   51.83   5.88     0       0     0; 
    5.94   5.88     10.28   0       0     0; 
    0        0          0       2.93   0     0;
    0        0          0       0       2.92 0;
    0        0          0       0       0   3.81]*1e9;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% properites of material layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'plate thickness [m]';                         lh = [1,1,1,1,1,1,1,1]*3.9/8/1000; % 3.9 mm total thickness
'stack angle sequence [deg]';              lalpha = [0 0 0 0 0 0 0 0];
'number of material layers';                 lay = length(lalpha);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single lamina mechanical properties (C tensor)
q11 = C(1,1); % [Pa]
q12 = C(1,2); % [Pa]
q22 = C(2,2); % [Pa]
%q13 = C(1,3); % [Pa] not present in the plate model
q44 = C(4,4); % [Pa]
q55 = C(5,5); % [Pa]
q66 = C(6,6); % [Pa]
rho = 1522.4; % [kg/m^3]
rho = repmat(rho,[lay,1]);

% damping coefficient proportional to mass matrix
etad_xy=0; % damping ratio
etad_z=0; % damping ratio
% damping coefficient at delamination region
etad_delam_xy=0; %
etad_delam_z=0; % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geometry definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V0=10; %[V] peak voltage 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% boundary conditions
% global b.c. vector with blocked degrees of freedom
BC=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pzt properties
% Noliac NCE51
 YE11=5.9e10; %[N/m^2] youngs modulus at constant electric field
 YE22=5.9e10; %[N/m^2] youngs modulus at constant electric field
 YE33=4.8e10; %[N/m^2] youngs modulus at constant electric field
 nu11=0.32; % Poisson ratio
% elastic compliance matrix
 Spzt=[ 1/YE11      -nu11/YE11    -nu11/YE33        0                      0                0;
       -nu11/YE11       1/YE11       -nu11/YE33        0                      0                0;
       -nu11/YE33    -nu11/YE33     1/YE33             0                      0                0;
           0                   0                  0         (1+nu11)/YE33           0                0;
           0                   0                  0                0               (1+nu11)/YE33     0;
           0                   0                  0                0                       0         (1+nu11)/YE11];%[m^2/N] efunda

% matrix of piezoelectric coupling constants - charge constants
      %sx    sy  sz  sxz syz sxy
dp=  [ 0     0       0     0   669 0;
         0     0       0    669  0   0;
      -208 -208   443   0     0   0]*10^-12;%[C/N]
% voltage constants
      %sx    sy  sz  sxz syz sxy
gp=  [ 0      0      0     0    38.9 0;
         0      0      0    38.9  0   0;
      -12.4 -12.4  26.3  0     0   0]*10^-3;%[Vm/N] 

% permittivity matrix (strain-charge form)
epsT=   [dp(1,5)/gp(1,5)       0                         0;
               0                 dp(1,5)/gp(1,5)            0;
               0                       0              dp(3,3)/gp(3,3)]; %[F/m]
% density
rho_pzt=7850;%[kg/m3]
pzt_thickness =0.5/1000; % pzt thickness [m]
theta_pzt = 0; % rotation angle of pzt [deg]