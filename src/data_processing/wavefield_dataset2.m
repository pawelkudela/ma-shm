clear all;close all;   warning off;clc;

load project_paths projectroot src_path;
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
%image_label_path = prepare_model_paths('raw','num','flat_shell','automesh_delam_rand'); % mesh parameters and labels
% prepare output paths
dataset_output_path = prepare_data_processing_paths('processed','num',modelname);
%figure_output_path = prepare_figure_paths(modelfolder,modelname);
%% Input for flat_shell
% load mesh parameters
%load([image_label_path,filesep,'mesh_parameters']);
%% input for post-processing
Nx=500;
Ny=500;
%shell_surface = {'top','bottom'}; % options: shell_surface = 'top'; shell_surface = 'bottom';
shell_surface = {'top'}; % options: shell_surface = 'top'; shell_surface = 'bottom';
field_variable='velocity';
motion=[3];%
%    motion - integer defining motion type: 
%    1) Ux
%    2) Uy
%    3) Uz
%    4) Fix
%    5) Fiy
%    6) Ux+h/2*UFix
%    7) Uy+h/2*UFiy
%    8) sqrt((Ux+h/2.*UFix).^2+(Uy+h/2.*UFiy).^2)
%    9) sqrt((Ux+h/2.*UFix).^2+(Uy+h/2.*UFiy).^2 + Uz.^2)
%tasks=[1:475];
tasks=[6];
for m=1%:2 % flat_shell_rand_1, flat_shell_rand_2 
% prepare input output paths
model_input_path = prepare_model_paths('interim','num','flat_shell',['flat_shell_rand_',num2str(m)]);      
        %%
    for test_case=tasks
        fprintf([modelname,' test case: %d\n'], test_case);
        for n=1:length(motion)
            [variable_name] = flat_shell_variable_names(field_variable,motion(n));
            for s=1:length(shell_surface)
                input_name = [model_input_path,filesep,num2str(test_case),'_output',filesep,'flat_shell_',variable_name,'_',num2str(test_case),'_',num2str(Nx),'x',num2str(Ny),shell_surface{s},'.mat'];
                processed_output_name = [dataset_output_path,filesep,num2str(test_case),'_output',filesep];
                %figure_output_name = [figure_output_path,num2str(test_case),'_output',filesep];
                filename=[processed_output_name,'1_flat_shell_',variable_name,'_',num2str(test_case),'_',num2str(Nx),'x',num2str(Ny),shell_surface{s}];
                if(overwrite||(~overwrite && ~exist([filename,'.png'],'file')))   
                    try
                        if ~exist(processed_output_name, 'dir')
                            mkdir(processed_output_name);
                        end
                        if exist( input_name , 'file')      
                        %% RUN DATA PROCESSING
                        load(input_name);% Data
                        [nx,ny,nft]=size(Data);
                        for k=1:nft
                            frame_filename=[processed_output_name,num2str(k),'_flat_shell_',variable_name,'_',num2str(test_case),'_',num2str(Nx),'x',num2str(Ny),shell_surface{s}];
                            frame=squeeze(Data(:,:,k));
                            F=frame2image(frame, frame_filename);
                        end
                        end
                    catch
                        fprintf('Failed test case no: %d\n', test_case);
                    end
                else
                    fprintf([modelname,' test case: %d %s %s already exist\n'], test_case,shell_surface{s},variable_name);
                end
            end
        end
    end
end



