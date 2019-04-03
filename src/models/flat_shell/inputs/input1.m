disp('.. Reading input data');
%% Signal definition
nft=1024*64;      % total number of samples
tt=0.5e-3;          % total calculation time [s] % 
t_1=0e-4;           % excitation initiation time [s]
f_1=50e3/5;        % frequency of the modulation signal [Hz]
f_2=5*f_1;          % frequency of the carrier signal [Hz]
nFrames=512;        % number of frames for animation
field_variable = 'velocity'; % field_variable for saving output data, string: 'displacement', 'velocity', 'acceleration' or 'all'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% properties of composite material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix(1) - aluminium
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Young modulus [Pa]';                      i_em(1) = 71e9;
'Poisson ratio';                               i_nim(1) = 0.33;
'density [kg/m3]';                           i_rhom(1) = 2700;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix(2) - A36 steel http://en.wikipedia.org/wiki/A36_steel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Young modulus [Pa]';                      i_em(2) = 200e9;
'Poisson ratio';                                i_nim(2) = 0.26;
'density [kg/m3]';                            i_rhom(2) = 7800;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix(3) - epoxy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Young modulus [Pa]';                      i_em(3) = 3.45e9;
'Poisson ratio';                                i_nim(3) = 0.35;
'density [kg/m3]';                            i_rhom(3) = 1097;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix(4) - pzt bonding layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Young modulus [Pa]';                      i_em(4) = 2e9;
'Poisson ratio';                               i_nim(4) = 0.34;
'density [kg/m3]';                           i_rhom(4) = 1170;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fibres(1) - fiber glass
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Young modulus [Pa]';                      i_ef(1) = 85e9; 
'Poisson ratio';                               i_nif(1) = 0.23; 
'density [kg/m3]';                          i_rhof(1) = 2150;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fibres(2) - fiber graphite/carbon
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'Young modulus [Pa]';                      i_ef(2) = 275.6e9; 
'Poisson ratio';                               i_nif(2) = 0.2; 
'density [kg/m3]';                          i_rhof(2) = 1900;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% properites of material layers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'plate thickness [m]';                         lh = [1,1,1,1]/1000;
'stack angle sequence [deg]';              lalpha = [0,90,90,0];
'stack matrix sequnce';                       lmat = [3,3,3,3];
'stack fibres sequece';                        lfib = [1,1,1,1];
'volume fraction of fibres';                   lvol=[0.192,0.192,0.192,0.192]; 
'number of material layers';                 lay = length(lalpha);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% damping coefficient proportional to mass matrix
etad_xy=0; % damping ratio
etad_z=0; % damping ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% geometry definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V0=10; %[V] peak voltage 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% boundary conditions
% global b.c. vector with blocked degrees of freedom
BC=[];
%% mesh
%meshfile='mesh/plate_Tomek_dens2_3mm1lay_pzt_mesh_2D.mat'; % 
% meshfile contains two matrices: coords and nodes and list of vectors with
% indices in local and global level
% element type: number of nodes in x, y and z direction respectively
% must be compatible with the mesh provided !
nx=6;ny=6;nz=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cd ..;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% delamination
% can be defined inside meshfile
delamination_layer = [2]; % delamination layer counting from the top after which the nodes are splitted (another column for second delamination)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% additional outputs
outputs=[25];

