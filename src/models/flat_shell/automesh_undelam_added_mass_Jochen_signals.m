clear all; close all;
% generate mesh for single geometry file (1 delamination, multi pzt)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Prepare output directories
% allow overwriting existing results if true
%overwrite=false;
overwrite=true;
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
% mesh parameters
shape_order = 5; % element shape function order, Number of nodes in one direction is shape_order+1
nPZT = 12; % number of piezoelectric trans
isAddedMassOn = true; % add nass instead of delamination
%%

mesh_filename = 'delam_Jochen_signals_D5_a_5mm_b_5mm_angle_0'; % added mass instead of delamination

figfilename = [figure_output_path,mesh_filename];

if(overwrite||(~overwrite && ~exist([mesh,filesep,mesh_filename,'added_mass.mat'], 'file')))
         %% RUN AUTOMESH
     try
        disp(mesh_filename);
     
        %[nodes,coords,den_under,den_above,I_G,I_L] = automesh_multi_pzt_simple(nPZT,shape_order,mesh_filename,modelfolder,figfilename,isDelamOn);
        [nodes,coords,addedMassEl,I_G,I_L] = automesh_multi_pzt_simple_added_mass(nPZT,shape_order,mesh_filename,modelfolder,figfilename,isAddedMassOn);

     catch
        fprintf(['Meshing failed:', mesh_filename,' \n']);
     end
else
    fprintf(['Mesh:', mesh_filename,' already exist\n']);
end
                
 