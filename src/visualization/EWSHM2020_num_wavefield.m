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
modelfolder = 'flat_shell'; % name of folder
modelname =  'flat_shell_EWSHM2020';
% prepare model output path

model_interim_path = prepare_model_paths('interim','num',modelfolder,modelname);
figure_output_path = prepare_figure_paths(modelfolder,modelname);
%% Input for flat_shell

tasks=[1:3];

%% input for post-processing

shell_surface = 'bottom'; % options: shell_surface = 'top'; shell_surface = 'bottom';
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
for test_case = tasks
interim_output_name = [model_interim_path,filesep,num2str(test_case),'_output',filesep];
figure_output_folder = [figure_output_path,num2str(test_case),'_output',filesep];

      if(overwrite||(~overwrite && ~exist([figure_output_folder], 'dir')))
            mkdir(figure_output_folder);
            load([interim_output_name,'flat_shell_Vz_1_982x982bottom']);
            Data_quarter=Data(1:491,1:491,:);% Data
            load([interim_output_name,'t_frames']); % t_frames (time vector)
            plot_meshgrid_frames_time(Data_quarter,t_frames,'bottom',test_case,selected_frames,figure_output_name,normalization,caxis_cut,ColorMapName,'velocity',3,fig_width,fig_height);

      else
          fprintf([modelname,' test case: %d already exist\n'], test_case);
      end

end



