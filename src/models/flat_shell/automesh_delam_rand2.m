clear all; close all;
% generate meshes for various delamination placements, sizes and orientations
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare output directories
% allow overwriting existing results if true
overwrite=false;
%overwrite=true;
% retrieve model name based on running file and folder
currentFile = mfilename('fullpath');
[pathstr,name,ext] = fileparts( currentFile );
idx = strfind( pathstr,filesep );
modelfolder = pathstr(idx(end)+1:end); % name of folder
modelname = name; 
% prepare model output path
image_label_path = prepare_model_paths('raw','num',modelfolder,modelname);
figure_output_path = prepare_figure_paths(modelfolder,modelname);
%%  INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=500; % image size in pixels
n=10; % mesh size of delamination positions
% input for automesh_delam
L=0.5; % plate length % for now fixed
W=0.5; % plate width % for now fixed

% range of delamination semi-major axis
amin=0.005;
amax=0.02;
% range of delamination semi-minor axis
bmin=0.005;
bmax=0.02;
% range of rotation angles
rotAngle_min = 0; 
rotAngle_max = 180;
r=0.005; % pzt radius
xpzt=0.25; % pzt x coordinate
ypzt=0.25; % pzt y coordinate

% mesh parameters
shape_order = 5; % element shape function order, Number of nodes in one direction is shape_order+1
CharacteristicLengthFactor = 0.068; 
CharacteristicLengthMin = 0.001; 
CharacteristicLengthMax = 0.2;
Smoothing = 1;
%%
% delamination scenario grid
[xn,yn]=meshgrid(0:n/(n-1):n,0:n/(n-1):n);
% x=xn*N/n; % coordinates in pixels
% y=yn*N/n; % coordinates in pixels
xc=xn*L/n; % coordinates in meters
yc=yn*W/n; % coordinates in meters
% plot(x,y,'ro');
% figure;
% plot(xc,yc,'ro');
mesh_parameters = struct('position_no',{},'xCenter',{},'yCenter',{},'a',{},'b',{},'rotAngle',{},'meshfile',{});
NofMeshes=476; %
% coordinates of delamination position
% I quarter
delta=0.01; % distance from pzt
a1=L/2+delta;
b1=W/2+delta;
xCenter1 = a1 + (L-a1).*rand(NofMeshes/4,1);
yCenter1 = b1 + (W-b1).*rand(NofMeshes/4,1);
% II quarter
a1=L/2+0.01;
b1=0;
xCenter2 = a1 + (L-a1).*rand(NofMeshes/4,1);
yCenter2 = b1 + (W/2-b1-delta).*rand(NofMeshes/4,1);
% III quarter
a1=0;
b1=0;
xCenter3 = a1 + (L/2-a1-delta).*rand(NofMeshes/4,1);
yCenter3 = b1 + (W/2-b1-delta).*rand(NofMeshes/4,1);
% IV quarter
a1=0;
b1=W/2+0.01;
xCenter4 = a1 + (L/2-a1-delta).*rand(NofMeshes/4,1);
yCenter4 = b1 + (W-b1).*rand(NofMeshes/4,1);
xCenter=[xCenter1;xCenter2;xCenter3;xCenter4];
yCenter=[yCenter1;yCenter2;yCenter3;yCenter4];

% mionor axis
b = bmin + (bmax-bmin).*rand(NofMeshes,1);
% major axis
a = amin + (amax-amin).*rand(NofMeshes,1);
% rotation angle
rotAngle = rotAngle_min + (rotAngle_max - rotAngle_min ).*rand(NofMeshes,1);
counter = 0;
% delamination position: all 4 quarters possible
for k=1:NofMeshes
    counter=counter+1;
    mesh_filename = ['m1_rand_single_delam_',num2str(counter)]; 
    mesh_parameters(counter).position_no = k;
    mesh_parameters(counter).xCenter = xCenter(k);
    mesh_parameters(counter).yCenter = yCenter(k);
    mesh_parameters(counter).a = a(k);
    mesh_parameters(counter).b = b(k);
    mesh_parameters(counter).rotAngle = rotAngle(k);
    mesh_parameters(counter).meshfile = mesh_filename;
end
for k=[257,318]
    mesh_filename = ['m1_rand_single_delam_',num2str(k)]; 
    figfilename = [figure_output_path,mesh_filename];
    image_label_filename = [image_label_path,filesep,mesh_filename];
    if(a(k) < b(k))
        atemp = a(k); a(k) = b(k); b(k) = atemp;
    end
    delam_image_label(N,xCenter(k)*N/L,yCenter(k)*N/W,a(k)*N/L,b(k)*N/W,rotAngle(k),image_label_filename);
    if(overwrite||(~overwrite && ~exist(['mesh',filesep,mesh_filename,'.mat'], 'file')))
             %% RUN AUTOMESH
         try
            disp(mesh_filename);
            [nodes,coords,den_under,den_above,I_G,I_L] = automesh_delam...
            (L,W,a(k),b(k),xCenter(k),yCenter(k),rotAngle(k),r,xpzt,ypzt,shape_order,CharacteristicLengthFactor,CharacteristicLengthMin,CharacteristicLengthMax,Smoothing,mesh_filename,modelfolder,figfilename);
         catch
            fprintf(['Meshing failed:', mesh_filename,' \n']);
         end
    else
        fprintf(['Mesh:', mesh_filename,' already exist\n']);
    end
end

save([image_label_path,filesep,'mesh_parameters'],'mesh_parameters');