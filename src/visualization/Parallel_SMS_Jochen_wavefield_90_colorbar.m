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
% prepare figure output path
modelname = name; 
figure_output_path = prepare_figure_paths(modelname);
modelfolder = 'flat_shell'; % name of folder
modelname =  'flat_shell_Jochen_wavefield_90';
% prepare model output path
model_interim_path = prepare_model_paths('interim','num',modelfolder,modelname);

%% Input for flat_shell

tasks=[1];

%% input for post-processing

shell_surface = 'top'; % options: shell_surface = 'top'; shell_surface = 'bottom';
%% input for figures
selected_frames = [64,72,100,108,120,128,256,512];
%selected_frames = [100,120,128];

%ColorMapName = 'jet';
ColorMapName = 'default';
%ColorMapName = 'map_sunset_interp'; % zero-symmetric
%ColorMapName = 'map_burd_interp'; % zero-symmetric
%ColorMapName = 'map_brewer_interp'; % orange-ish
%ColorMapName = 'map_iridescent_interp'; % blue-ish
%ColorMapName = 'map_discrete_rainbow'; % not suitable for wavefields
caxis_cut = 0.8;
normalization = false; % normalization to the highest value of the wavefield
fig_width =5; % figure widht in cm
fig_height=5+2; % figure height in cm (+2 for colorbar)
%    ColorMapName - custom map name: 'map_sunset', 'map_sunset_interp', 'map_burd', 'map_burd_interp'
%    'map_brewer', 'map_brewer_interp', 'map_iridescent', 'map_iridescent_interp', 'map_discrete_rainbow'
%%
test_case =1;
scaling=20;
interim_output_name = [model_interim_path,filesep,num2str(test_case),'_output',filesep];
figure_output_name = [figure_output_path,num2str(test_case),'_output',filesep];
 if ~exist(figure_output_name, 'dir')
    mkdir(figure_output_name);
 end
load([interim_output_name,'flat_shell_Vz_1_500x500top']);
Data_quarter=scaling*Data(11:250,11:250,:);
plot_meshgrid_frames_colorbar(Data_quarter,shell_surface,test_case,selected_frames,[figure_output_name,'colorbar_'],normalization,caxis_cut,ColorMapName,'velocity',3,fig_width,fig_height);
% caxis cut
caxis_cut = 0.1;
Data_zoom=scaling*Data(11:250,11:250,:);
plot_meshgrid_frames_colorbar(Data_zoom,'top',test_case,selected_frames,[figure_output_name,'colorbar_cut_'],normalization,caxis_cut,ColorMapName,'velocity',3,fig_width,fig_height);
caxis_cut = 0.8;
load([interim_output_name,'flat_shell_Vz_1_500x500bottom']);
Data_quarter=scaling*Data(11:250,11:250,:);
plot_meshgrid_frames_colorbar(Data_quarter,'bottom',test_case,selected_frames,[figure_output_name,'colorbar_'],normalization,caxis_cut,ColorMapName,'velocity',3,fig_width,fig_height);


