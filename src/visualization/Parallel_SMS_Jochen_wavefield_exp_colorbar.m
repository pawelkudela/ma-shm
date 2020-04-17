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

%% input for experiment
list = {'240px240p_50kHz'};  
% prepare exp input path
input_data_path = fullfile( projectroot, 'data','interim','exp', filesep );

%% input for post-processing

shell_surface = 'bottom'; % options: shell_surface = 'top'; shell_surface = 'bottom';
%% input for figures
selected_frames = [257,401];


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
scaling=1;

figfilename = list{test_case};
load([input_data_path,figfilename]);
Data=permute(Data,[2,1,3]);
plot_meshgrid_frames_colorbar(Data,shell_surface,test_case,selected_frames,[figure_output_path,'colorbar_'],normalization,caxis_cut,ColorMapName,'velocity',3,fig_width,fig_height);
