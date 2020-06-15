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
modelfolder = pathstr(idx(end)+1:end); % name of folder
modelname = name; 
%image_label_path = prepare_model_paths('raw','num','flat_shell','automesh_delam_rand'); % mesh parameters and labels
image_label_path=fullfile(filesep,filesep, 'odroid-sensors','sensors','aidd','data','raw','num','dataset2_labels_out',filesep);
% prepare output paths
dataset_output_path = prepare_data_processing_paths('processed','num',modelname);
figure_output_path = prepare_figure_paths(modelname);
%% input for figures
Cmap = jet(256); 
caxis_cut = 0.8;
fig_width =5; % figure widht in cm
fig_height=5; % figure height in cm
%% Input for signal processing

WL = [0.5 0.5];
mask_thr = 0.85;
PLT = 0.5;
%threshold = 0.012; % threshold for binarization
%threshold = 0.008; % threshold for binarization
threshold = 0.0075; % threshold for binarization
%% input for post-processing
Nx=500;
Ny=500;
shell_surface = {'top','bottom'}; % options: shell_surface = 'top'; shell_surface = 'bottom';
%shell_surface = {'bottom'}; 
%shell_surface = {'top'}; 
field_variable='velocity';
%motion=[3,8];
motion=[3];
tasks=[1:475];
%tasks=[16];
IoU=zeros(length(tasks),length(shell_surface),length(motion));
area=zeros(length(tasks),length(shell_surface),length(motion));
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
                if(s == 1) threshold = 0.012; end
                if(s == 2) threshold = 0.008; end
                processed_output_name = [dataset_output_path,filesep];
                %figure_output_name = [figure_output_path,num2str(test_case),'_output',filesep];
                fname=['ERMSF_flat_shell_',variable_name,'_',num2str(test_case),'_',num2str(Nx),'x',num2str(Ny),shell_surface{s}];
                filename=[processed_output_name,fname];
                if(overwrite||(~overwrite && ~exist([filename,'.mat'],'file')))   
                    try
                        if ~exist(processed_output_name, 'dir')
                            mkdir(processed_output_name);
                        end
                        if exist( input_name , 'file')    
                            load([model_input_path,filesep,num2str(test_case),'_output',filesep,'t_frames']);
                            time=t_frames(1:512);
                        %% RUN DATA PROCESSING
                        load(input_name);% Data
                        [nx,ny,nft]=size(Data);
                        
                        %% Adaptive filtering
                        [RMSF,ERMSF,WRMSF] = AdaptiveFiltering(Data,time,WL,mask_thr,PLT);
                       %% save picture
                       figure;
                       imagesc(ERMSF);
                       colormap(Cmap);
                       set(gca,'YDir','normal');axis square;axis off;
                       set(gcf,'color','white');
                       %Smin=min(min(ERMSF));
                       Smin=0;
                       Smax=max(max(ERMSF));
                       set(gcf,'Renderer','zbuffer');
                       caxis([caxis_cut*Smin,caxis_cut*Smax]);
                        set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
                        set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
                        % remove unnecessary white space
                        %set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
                        set(gcf,'PaperPositionMode','auto');

                        print([figure_output_path,fname],'-dpng', '-r600'); 
                       %% binarization                
                        Binary = uint8(ERMSF >= threshold);
                        save(filename,'ERMSF','WL','mask_thr','Binary','threshold');
                        % plot
                        whitethreshold = .05;
                        blackthreshold = .05;
                        CmapB = 1-([blackthreshold:1/255:1-whitethreshold ; blackthreshold:1/255:1-whitethreshold ; blackthreshold:1/255:1-whitethreshold]');
                        figure
                        imagesc(Binary)
                        colormap(1-CmapB)
                        set(gca,'YDir','normal');
                        axis square;axis off;
                        set(gcf,'color','white');
                        set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
                        set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
                        print([figure_output_path,'Binary_',fname],'-dpng', '-r600'); 
                        %% intersection over union
                        labelname=[image_label_path,'m1_rand_single_delam_',num2str(test_case),'.png'];
                        A=imread(labelname)/255;
                        IoU(test_case,s,n)=intersect_over_union_fun(flipud(Binary),A);
                        area(test_case,s,n)=sum(sum(Binary))/(Nx*Ny)*WL(1)*1e3*WL(2)*1e3; % [mm^2]
                        disp('Intersection over union: ');IoU(test_case,s,n)
                        close all;
                        % area ?
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
save([processed_output_name,'IoU'],'IoU','area');


