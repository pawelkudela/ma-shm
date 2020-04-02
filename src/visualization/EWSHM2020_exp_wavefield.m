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
% create path to the experimental data folder
input_data_path = fullfile( projectroot, 'data','raw','exp', filesep );
% prepare figure output path
figure_output_path = prepare_figure_paths(modelname);

%% input for figures
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
% full wavefield measurements
list = {'491x491p_16_5kHz_5HC_x5_15Vpp', ...     
         '492x492p_50kHz_5HC_x20_15Vpp', ...
         '491x491p_100kHz_5HC_x20_15Vpp'};    
frame_list ={[16,32,64:64:512],[16,32,64:64:512]*2,[16,32,64:64:512]*2};
%%
tasks =1:3; % file numbers from list
folder  = input_data_path;
nFile   = length(tasks);
for test_case =tasks
    fprintf([modelname,' test case: %d\n'], test_case);
    filename = list{test_case};
    selected_frames = frame_list{test_case};
    figure_output_folder = [figure_output_path,num2str(test_case),'_output',filesep];
     figfilename = [figure_output_folder,filename,'_','frame',num2str(selected_frames(1))];
      if(overwrite||(~overwrite && ~exist([figure_output_folder], 'dir')))
            mkdir(figure_output_folder);
            load([input_data_path,filename]); % Data, time
            plot_exp_frames(Data,time,filename,selected_frames,figure_output_folder,normalization,caxis_cut,ColorMapName,fig_width,fig_height);
      else
          fprintf([modelname,' test case: %d already exist\n'], test_case);
      end
end




