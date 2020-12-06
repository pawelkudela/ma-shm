disp('.. Reading input data');
%% Signal definition
tt= 0.00199609375;         % total calculation time [s] % tt=time(512)
nFrames_exp = 512; % number of time steps taken from experiment
dt_exp = tt/(nFrames_exp-1); % time step in experiment
division = 2^8; % number of divisions of experimental time step
dt_num = dt_exp/division; % time step in numerical simulation
nft = (nFrames_exp-1)*division+1; % total number of samples
t_1=0e-4;           % excitation initiation time [s]
f_1=16.5e3/5;        % frequency of the modulation signal [Hz]
f_2=5*f_1;          % frequency of the carrier signal [Hz]
nFrames=512;        % number of frames for animation
frame_no = ((([1:nFrames]-1)*division)*nFrames_exp/nFrames+1);
field_variable = 'velocity'; % field_variable for saving output data, string: 'displacement', 'velocity', 'acceleration' or 'all'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% properties of composite material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix(1) - aluminium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Young modulus [Pa]';                      i_em(1) = 71e9;
% 'Poisson ratio';                               i_nim(1) = 0.33;
% 'density [kg/m3]';                           i_rhom(1) = 2700;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % matrix(2) - A36 steel http://en.wikipedia.org/wiki/A36_steel
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Young modulus [Pa]';                      i_em(2) = 200e9;
% 'Poisson ratio';                                i_nim(2) = 0.26;
% 'density [kg/m3]';                            i_rhom(2) = 7800;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % matrix(3) - epoxy
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Young modulus [Pa]';                      i_em(3) = 3.45e9;
% 'Poisson ratio';                                i_nim(3) = 0.35;
% 'density [kg/m3]';                            i_rhom(3) = 1097;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % matrix(4) - pzt bonding layer
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Young modulus [Pa]';                      i_em(4) = 2e9;
% 'Poisson ratio';                               i_nim(4) = 0.34;
% 'density [kg/m3]';                           i_rhom(4) = 1170;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % matrix(5) - semi-optimized properties
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Young modulus [Pa]';                      i_em(5) = 3.43e9;
% 'Poisson ratio';                               i_nim(5) = 0.35;
% 'density [kg/m3]';                           i_rhom(5) = 1250;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fibres(1) - fiber glass
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Young modulus [Pa]';                      i_ef(1) = 85e9; 
% 'Poisson ratio';                               i_nif(1) = 0.23; 
% 'density [kg/m3]';                          i_rhof(1) = 2150;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fibres(2) - fiber graphite/carbon
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Young modulus [Pa]';                      i_ef(2) = 275.6e9; 
% 'Poisson ratio';                               i_nif(2) = 0.2; 
% 'density [kg/m3]';                          i_rhof(2) = 1900;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % fibres(3) - carbon fiber Hexply
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 'Young modulus [Pa]';                      i_ef(3) = 211.2e9; 
% 'Poisson ratio';                               i_nif(3) = 0.2; 
% 'density [kg/m3]';                          i_rhof(3) = 1900;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% properites of material layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'plate thickness [m]';                         lh = [1,1,1,1,1,1,1,1]*2.85/8/1000; % 2.85 mm total thickness
'stack angle sequence [deg]';              lalpha = [90,90,90,90,90,90,90,90];
% 'stack matrix sequnce';                       lmat = [5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5];
% 'stack fibres sequece';                        lfib = [3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3];
% 'volume fraction of fibres';                   lvol=[0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44]; 
'number of material layers';                 lay = length(lalpha);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% single lamina mechanical properties (C tensor)
% q11 = 139.3e9; % [Pa]
% q12 = 6.73e9; % [Pa]
% q22 = 12.64e9; % [Pa]
% q44 = 3.07e9; % [Pa]
% q55 = 5.15e9; % [Pa]
% q66 = 4.16e9; % [Pa]
% with applied C32_2_C2D
q11 = 135.72e9; % [Pa]
q12 = 3.35e9; % [Pa]
q22 = 9.44e9; % [Pa]
q44 = 3.07e9; % [Pa]
q55 = 5.15e9; % [Pa]
q66 = 4.16e9; % [Pa]
%
rho = 1574.1; % [kg/m^3]
rho = repmat(rho,[lay,1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% damping coefficient proportional to mass matrix
etad_xy=0; % damping ratio
etad_z=0; % damping ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geometry definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V0=100; %[V] peak voltage 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% boundary conditions
% global b.c. vector with blocked degrees of freedom
BC=[];
%% delamination
% can be defined inside meshfile
delamination_layer = [4]; % delamination layer counting from the top after which the nodes are splitted (another column for second delamination)
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