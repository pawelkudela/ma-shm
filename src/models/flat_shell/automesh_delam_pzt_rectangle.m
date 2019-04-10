clear all; close all;
% generate mesh for single geometry file (1 delamination, square pzt)
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

%%

mesh_filename = 'rectangle_pzt_test3'; 

figfilename = [figure_output_path,mesh_filename];

if(overwrite||(~overwrite && ~exist([mesh,filesep,mesh_filename,'.mat'], 'file')))
         %% RUN AUTOMESH
     try
        disp(mesh_filename);
        [nodes,coords,den_under,den_above,I_G,I_L] = automesh_simple(shape_order,mesh_filename,modelfolder,figfilename);
     catch
        fprintf(['Meshing failed:', mesh_filename,' \n']);
     end
else
    fprintf(['Mesh:', mesh_filename,' already exist\n']);
end
                
 