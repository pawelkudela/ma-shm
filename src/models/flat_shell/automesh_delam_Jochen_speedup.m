clear all; close all;
% generate mesh for single geometry file (1 delamination, multi pzt discs)
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
nPZT = 1; % number of piezoelectric trans
isDelamOn = true; % delamination case
%%

mesh_filename = {'delam_Jochen_wavefield_a_5mm_b_5mm_angle_0_L_300_W_300';
                        'delam_Jochen_wavefield_a_5mm_b_5mm_angle_0_L_400_W_400';
                        'delam_Jochen_wavefield_a_5mm_b_5mm_angle_0_L_500_W_500';
                        'delam_Jochen_wavefield_a_5mm_b_5mm_angle_0_L_700_W_700';
                        'delam_Jochen_wavefield_a_5mm_b_5mm_angle_0_L_900_W_900';
                        'delam_Jochen_wavefield_a_5mm_b_5mm_angle_0_L_1000_W_1000'}; 
                    
for k=1:length(mesh_filename)
    figfilename = [figure_output_path,mesh_filename{k}];

    if(overwrite||(~overwrite && ~exist([mesh,filesep,mesh_filename{k},'.mat'], 'file')))
             %% RUN AUTOMESH
         try
            disp(mesh_filename);

            [nodes,coords,den_under,den_above,I_G,I_L] = automesh_multi_pzt_simple(nPZT,shape_order,mesh_filename{k},modelfolder,figfilename,isDelamOn);
         catch
            fprintf(['Meshing failed:', mesh_filename,' \n']);
         end
    else
        fprintf(['Mesh:', mesh_filename,' already exist\n']);
    end
end                
 