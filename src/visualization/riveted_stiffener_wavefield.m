clear all;close all;   warning off;clc;

load project_paths projectroot src_path;
%% Prepare output directories
% allow overwriting existing results if true
overwrite=false;
overwrite=true;
% retrieve model name based on running file and folder
currentFile = mfilename('fullpath');
[pathstr,name,ext] = fileparts( currentFile );
idx = strfind( pathstr,filesep );
modelfolder = 'solid3D'; % name of folder
modelname =  'main_spec3Dcuda_global_cpu_gpu3_wave_sens_opt_single_PCT_rivet';
modelname_out =  name;
% prepare model output path

model_interim_path = prepare_model_paths('interim','num',modelfolder,modelname);
figure_output_path = prepare_figure_paths(modelfolder,modelname_out);
%% Input for flat_shell

tasks=[7:8];

%% input for post-processing

shell_surface = 'bottom'; % options: shell_surface = 'top'; shell_surface = 'bottom';
%% input for figures
selected_frames = [64,128,160,192,250];
%ColorMapName = 'jet';
ColorMapName = 'default';
%ColorMapName = 'map_sunset_interp'; % zero-symmetric
%ColorMapName = 'map_burd_interp'; % zero-symmetric
%ColorMapName = 'map_brewer_interp'; % orange-ish
%ColorMapName = 'map_iridescent_interp'; % blue-ish
%ColorMapName = 'map_discrete_rainbow'; % not suitable for wavefields

normalization = false; % normalization to the highest value of the wavefield
fig_width =5; % figure widht in cm
fig_height=5; % figure height in cm
%    ColorMapName - custom map name: 'map_sunset', 'map_sunset_interp', 'map_burd', 'map_burd_interp'
%    'map_brewer', 'map_brewer_interp', 'map_iridescent', 'map_iridescent_interp', 'map_discrete_rainbow'
%%
for test_case = tasks
    ColorMapName = 'default';
    caxis_cut = 0.2;
interim_output_name = [model_interim_path,filesep,num2str(test_case),'_output',filesep];
figure_output_folder = [figure_output_path,num2str(test_case),'_output',filesep];

      if(overwrite||(~overwrite && ~exist([figure_output_folder], 'dir')))
            mkdir(figure_output_folder);
            load([interim_output_name,'solid3D_Uz_',num2str(test_case),'_500x500bottom']);
            load([interim_output_name,'t_frames']); % t_frames (time vector)
            plot_meshgrid_frames_time_3D(Data,t_frames,'bottom',test_case,selected_frames,figure_output_folder,normalization,caxis_cut,ColorMapName,'displacement',3,fig_width,fig_height);
            
            ColorMapName = 'jet';caxis_cut = 0.5;
            load([interim_output_name,'solid3D_in_plane_',num2str(test_case),'_500x500bottom']);      
            plot_meshgrid_frames_time_3D(Data,t_frames,'bottom',test_case,selected_frames,figure_output_folder,normalization,caxis_cut,ColorMapName,'displacement',4,fig_width,fig_height);

      else
          fprintf([modelname,' test case: %d already exist\n'], test_case);
      end

end



