clear all;close all;   warning off;clc;

load project_paths projectroot src_path;
%% Prepare output directories
% allow overwriting existing results if true
%overwrite=false;
overwrite=true;
% figure parameters
% size 12cm by 8cm (1-column text)
%fig_width = 12; fig_height = 8; 
% size 7cm by 5cm (2-column text)
fig_width = 7; fig_height = 7; 
modelfolder = 'flat_shell'; % name of folder
caxis_cut = 0.5;
%% Input for flat_shell
% load mesh parameters
%load([image_label_path,filesep,'mesh_parameters']);
%% input for post-processing
Nx=500;
Ny=500;
shell_surface = {'top'}; % options: shell_surface = 'top'; shell_surface = 'bottom';
field_variable='velocity';
motion=[3];
tasks=[1];
%tasks=[1];
j=0;
for m=[1,6,11,16] % output no
    % prepare output paths
    modelname = ['flat_shell_Jochen_random_exc_signals',num2str(m)];
    figure_output_path = prepare_figure_paths(modelfolder,modelname);
% prepare input output paths
    model_input_path = prepare_model_paths('interim','num','flat_shell',['flat_shell_Jochen_random_exc_signals',num2str(m)]);      
        %%
    for test_case=tasks
        j=j+1;
        fprintf([modelname,' test case: %d\n'], test_case);
        for n=1:length(motion)
            [variable_name] = flat_shell_variable_names(field_variable,motion(n));
            for s=1:length(shell_surface)
                input_name = [model_input_path,filesep,num2str(test_case),'_output',filesep,'flat_shell_',variable_name,'_',num2str(test_case),'_',num2str(Nx),'x',num2str(Ny),shell_surface{s},'.mat'];
                figure_output_name = [figure_output_path,num2str(test_case),'_output',filesep];
                filename=[figure_output_name,'WRMS_flat_shell_',variable_name,'_',num2str(test_case),'_',num2str(Nx),'x',num2str(Ny),shell_surface{s}];
                if(overwrite||(~overwrite && ~exist([filename,'.png'],'file')))   
                    try
                        if ~exist(figure_output_name, 'dir')
                            mkdir(figure_output_name);
                        end
                        if exist( input_name , 'file')      
                        %% RUN DATA PROCESSING
                        load(input_name);% Data
                        [nx,ny,nft]=size(Data);
                        nft=256;
                        %RMS = abs(sqrt(sum(Data(2:end-1,2:end-1,1:nft).^2,3)));
                        Weighted_Data = zeros(nx,ny,nft);
                        for frame=1:nft
                             Data(:,:,frame) = mymedian3x3(Data(:,:,frame)); % 3x3 median filtering
                            Weighted_Data(:,:,frame) = Data(:,:,frame)*sqrt(frame^3);
                        end
                        WRMS = abs(sqrt(sum(Weighted_Data(2:end-1,2:end-1,1:end).^2,3)));
                        h=figure(j);
                        colormap jet;
                        imagesc(WRMS); 
                        caxis([0,caxis_cut*max(max(WRMS))]);
                        %imagesc(RMS); 
                        %caxis([0,caxis_cut*max(max(RMS))]);
                        set(h,'Color','w');
                        set(gca,'YDir','normal'); 
                        axis off;axis square;
                        
                        set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
                        set(h, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); % 
                        % remove unnecessary white space
                        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
                        h.PaperPositionMode   = 'auto';
                        print(filename,'-dpng', '-r600'); 
                        end
                        close all;
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



