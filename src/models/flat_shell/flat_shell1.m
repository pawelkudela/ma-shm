clear all;close all;   warning off;clc;

load project_paths projectroot src_path;
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
model_output_path = prepare_model_paths('raw','num',modelfolder,modelname);
model_interim_path = prepare_model_paths('interim','num',modelfolder,modelname);
%% Input for flat_shell
tasks=[6];
mode='gpu'; % options: mode='cpu';mode='gpu';
meshfile='mesh/plate_Tomek_dens2_3mm1lay_pzt_mesh_2D.mat'; % 
%% input for post-processing
Nx=500;
Ny=500;
for test_case=tasks
    output_name = [model_output_path,filesep,num2str(test_case),'_output',filesep];
    interim_output_name = [model_interim_path,filesep,num2str(test_case),'_output',filesep];
    if(overwrite||(~overwrite && ~exist(output_name, 'dir')))
        fprintf([modelname,' test case: %d\n'], test_case);
        try
            if ~exist(output_name, 'dir')
                mkdir(output_name);
            end
            if ~exist(interim_output_name, 'dir')
                mkdir(interim_output_name);
            end
            main_flat_shell(test_case,meshfile,mode,output_name,tasks);
            [Data] = spec2meshgrid_flat_shell(test_case,Nx,Ny,'velocity',3,'upper',interim_output_name); % Vz 
            [Data] = spec2meshgrid_flat_shell(test_case,Nx,Ny,'velocity',8,'upper',interim_output_name); % sqrt((Vx+h/2.*VFix).^2+(Vy+h/2.*VFiy).^2)
            delete([output_name,'*frame*']); % delete frames
        catch
            fprintf('Failed test case no: %d\n', test_case);
        end
    else
        fprintf([modelname,' test case: %d already exist\n'], test_case);
    end
end




