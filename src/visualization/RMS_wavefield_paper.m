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
image_label_path = prepare_model_paths('raw','num','flat_shell','automesh_delam_rand'); % mesh parameters and labels
% prepare output paths
dataset_output_path = prepare_data_processing_paths('processed','num',modelname);
% prepare output paths
figure_output_path = prepare_figure_paths(modelname);
%% input for figures
Cmap = gray(256); 
caxis_cut = 0.8;
fig_width =5; % figure widht in cm
fig_height=5; % figure height in cm
%% Input for flat_shell
% load mesh parameters
%load([image_label_path,filesep,'mesh_parameters']);
%% input for post-processing
Nx=500;
Ny=500;
shell_surface = {'top','bottom'}; % options: shell_surface = 'top'; shell_surface = 'bottom';
field_variable='velocity';
motion=[3,8];
tasks=[1,433,438,454];
%tasks=[1];
for m=1:3 % flat_shell_rand1, flat_shell_rand2 
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
                filename=[processed_output_name,'RMS_flat_shell_',variable_name,'_',num2str(test_case),'_',num2str(Nx),'x',num2str(Ny),shell_surface{s}];
                filename2=['RMS_flat_shell_',variable_name,'_',num2str(test_case),'_',num2str(Nx),'x',num2str(Ny),shell_surface{s}];
  
                if(overwrite||(~overwrite && ~exist([figure_output_path,filename2,'.png'],'file')))   
                    try
                        if ~exist(processed_output_name, 'dir')
                            mkdir(processed_output_name);
                        end
                        if exist( input_name , 'file')      
                        %% RUN DATA PROCESSING
                        load(input_name);% Data
                        [nx,ny,nft]=size(Data);
                        RMS = abs(sqrt(sum(Data(:,:,1:256).^2,3)));
                        A=rms2image(RMS, filename);
                %% save picture
                       figure;
                       imagesc(A);
                       colormap(Cmap);
                       %set(gca,'YDir','normal');
                       axis square;axis off;
                       set(gcf,'color','white');

                       Smin=0;
                       Smax=max(max(A));
                       set(gcf,'Renderer','zbuffer');
                       caxis([caxis_cut*Smin,caxis_cut*Smax]);
                        set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
                        set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
                        % remove unnecessary white space
                        %set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
                        set(gcf,'PaperPositionMode','auto');

                        print([figure_output_path,filename2],'-dpng', '-r600'); 
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


