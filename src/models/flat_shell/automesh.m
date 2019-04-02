% testing automeshing with gmsh



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=100; % image size in pixels
n=10; % mesh size of delamination positions
% input for automesh_delam
L=0.5; % plate length % for now fixed
W=0.5; % plate width % for now fixed
a=0.02; % delamination semi-major axis
b=0.01; % delamination semi-minor axis
xCenter=0.5;   % delamination x coordinate
yCenter=0.5;   % delamination y coordinate
rotAngle = 30; % delamination rotation [0:180)
r=0.005; % pzt radius
xpzt=0.25; % pzt x coordinate
ypzt=0.25; % pzt y coordinate

% mesh parameters
shape_order = 5; % element shape function order, Number of nodes in one direction is shape_order+1
CharacteristicLengthFactor = 0.08; 
CharacteristicLengthMin = 0.001; 
CharacteristicLengthMax = 0.2;
Smoothing = 1;
% filenames and paths
mesh_filename = 'delam_1_test';
modelfolder = 'flat_shell';
%%
% i=8;j=4; % delamination grid point (i -  rows  , j - columns)
% m=(j-1)*n+i; % example number
% % delamination scenario grid
% [xn,yn]=meshgrid(0:n/(n-1):n,0:n/(n-1):n);
% x=xn*N/n; % coordinates in pixels
% y=yn*N/n; % coordinates in pixels
% xc=xn*L/n; % coordinates in meters
% yc=yn*W/n; % coordinates in meters
% % coordinates of delaminatio position
% xCenter = xc(j,i);
% yCenter = yc(j,i);
% % plot(x,y,'ro');
% % figure;
% % plot(xc,yc,'ro');
% return;
%% RUN AUTOMESH
[nodes,coords,den_under,den_above,IG1,IG2,IG3,IG4,IG5,IG6,IG7,IG8,IG9,IG10,IG11,IG12,IL1,IL2,IL3,IL4,IL5,IL6,IL7,IL8,IL9,IL10,IL11,IL12] = automesh_delam...
    (L,W,a,b,xCenter,yCenter,rotAngle,r,xpzt,ypzt,shape_order,CharacteristicLengthFactor,CharacteristicLengthMin,CharacteristicLengthMax,Smoothing,mesh_filename,modelfolder);

%%



% save([spec_mesh_output_path,'pztEl_', mesh_filename,'.mat'],'pztEl');
% save([spec_mesh_output_path,'delamEl_', mesh_filename,'.mat'],'delamEl');
% % pause;
% close all;
% figfilename = 'test_mesh';
% print([figfilename],'-dpng', '-r600'); 