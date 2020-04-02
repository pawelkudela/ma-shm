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
figure_output_path = prepare_figure_paths(modelfolder,modelname);
%% Input for flat_shell

% input for constant parameters
input_file_no = [21,22,23]; % frequency = [16.5e3/5,50e3/5,100e3/5];
tasks=[1:3];
mode='gpu'; % options: mode='cpu';mode='gpu';
meshfile=fullfile('mesh','single_pzt_single_delam_large_plate_EWSHM2020');
%% input for post-processing
Nx=982;
Ny=982;
shell_surface = 'top'; % options: shell_surface = 'top'; shell_surface = 'bottom';
%% input for figures
selected_frames = [16,32,64:64:512];
ColorMapName = 'jet';
%ColorMapName = 'default';
%ColorMapName = 'map_sunset_interp'; % zero-symmetric
%ColorMapName = 'map_burd_interp'; % zero-symmetric
%ColorMapName = 'map_brewer_interp'; % orange-ish
%ColorMapName = 'map_iridescent_interp'; % blue-ish
%ColorMapName = 'map_discrete_rainbow'; % not suitable for wavefields
caxis_cut = 0.8;
normalization = false; % normalization to the highest value of the wavefield
fig_width =5; % figure widht in cm
fig_height=5; % figure height in cm
%    ColorMapName - custom map name: 'map_sunset', 'map_sunset_interp', 'map_burd', 'map_burd_interp'
%    'map_brewer', 'map_brewer_interp', 'map_iridescent', 'map_iridescent_interp', 'map_discrete_rainbow'
%%
for test_case=tasks
    input_no = input_file_no(test_case);
    output_name = [model_output_path,filesep,num2str(test_case),'_output',filesep];
    interim_output_name = [model_interim_path,filesep,num2str(test_case),'_output',filesep];
    figure_output_name = [figure_output_path,num2str(test_case),'_output',filesep];
    if(overwrite||(~overwrite && ~exist([output_name,'time',num2str(test_case),'.mat'],'file')))
        fprintf([modelname,' test case: %d\n'], test_case);
        try
            if ~exist(output_name, 'dir')
                mkdir(output_name);
            end
            if ~exist(interim_output_name, 'dir')
                mkdir(interim_output_name);
            end
             if ~exist(figure_output_name, 'dir')
                mkdir(figure_output_name);
             end
             
            %% RUN MODEL
            actuator_no = 1;
            t_frames=main_flat_shell_multi_pzt_c(actuator_no,test_case,input_no,meshfile,mode,output_name,tasks);
            %
            t_frames_filename=fullfile(interim_output_name,'t_frames');
            save(t_frames_filename,'t_frames');
            %% RUN POSTPROCESSING
            % out-of-plane
            
            [Data] = spec2meshgrid_flat_shell(test_case,input_no,meshfile,Nx,Ny,'velocity',3,'bottom',output_name,interim_output_name); % Vz 
            %plot_meshgrid_frames(Data(1:491,1:491,:),'bottom',test_case,selected_frames,figure_output_name,normalization,caxis_cut,ColorMapName,'velocity',3,fig_width,fig_height);
            plot_meshgrid_frames_time(Data(1:491,1:491,:),t_frames,'bottom',test_case,selected_frames,figure_output_name,normalization,caxis_cut,ColorMapName,'velocity',3,fig_width,fig_height);
            delete([output_name,'*frame*']); % delete frames
        catch
            fprintf('Failed test case no: %d\n', test_case);
        end
    else
        fprintf([modelname,' test case: %d already exist\n'], test_case);
    end
end




