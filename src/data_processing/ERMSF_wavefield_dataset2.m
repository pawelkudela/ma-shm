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
%figure_output_path = prepare_figure_paths(modelfolder,modelname);
%% Input for signal processing

Nmed = 3; k = 3; Vh = 1014; Vv = 923;  WL = [0.44 0.44]; timeframe = 1:130;

%NFrames = 1:30; thresh = 5;
%NFrames = 1:8; thresh = 8;
%NFrames = 1:5; thresh = 5;
%NFrames = 1:5; thresh = 8;
%NFrames = 1:3; thresh = 5;
%NFrames = 1:3; thresh = 10;
%NFrames = 1:5; thresh = 10;
%NFrames = 1:30; thresh = 10;
%NFrames = 1:20; thresh = 10;
%NFrames = 1:20; thresh = 8;
%NFrames = 1:20; thresh = 5;
%NFrames = 1:12; thresh = 10;
%NFrames = 1:12; thresh = 12;
%NFrames = 1:12; thresh = 8; % detect edges
%NFrames = 1:12; thresh = 6;
%NFrames = 1:12; thresh = 4; % edges are blurred
%NFrames = 1:10; thresh = 12;
%NFrames = 1:26; thresh = 6;
NFrames = 1:512; thresh = 10;
PRT = 1;% if PRT = 1 then print
%a=64; % no of points for extension beyond the plate dimensions (zero padding)
a=0;
%a=8;
%a=4;
%% input for post-processing
Nx=500;
Ny=500;
%shell_surface = {'top','bottom'}; % options: shell_surface = 'top'; shell_surface = 'bottom';
shell_surface = {'bottom'}; 
%shell_surface = {'top'}; 
field_variable='velocity';
%motion=[3,8];
motion=[3];
tasks=6;%[1:475];
%tasks=[1];

for m=1:1 % flat_shell_rand1, flat_shell_rand2 
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
                if(overwrite||(~overwrite && ~exist([filename,'.png'],'file')))   
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
                        % medial filtering not needed for numerical data
                        if (a>0)
                            FData=zeros(nx+a,ny+a,nft);
                            FData(a/2:nx+a/2-1,a/2:ny+a/2-1,:)=Data;
                        else
                            FData=Data;
                        end
                        %% EData (Amlidute composated data)
                        %trs = 0.8; %numerical 0.8 %experimentak 0.95
                        %[EData,E,P,E2] = EqualizeData(FData,time,Vv,Vh,trs);
                        %% Energy plot
                        %End_cut = 5;
                        %[hE,strt,endd] = propPlotEnery(time,E,E2,End_cut,trs,PRT,['Energy_c' num2str(c)]);
                       
                        %% FFT 2D
                        [FFTData] = propFFT2(FData);
                        %% Adaptive filter mask
                       %[FilterMask] = propMask(FFTData,NFrames,thresh);           
                       [FilterMask,CDF_x,CDF] = propMask_cum(FFTData,NFrames,thresh);
                      %% Gaussian blur for the mask
                      sigma = 30;
                      %sigma = 10;
                      %[M_G] = LowPassMaskGauss(FilterMask,sigma); % Gaussian blur
                      %FilterMask = FilterMask.*(1-M_G);
                      %FilterMask = filter2( Gauss(5,1),FilterMask);
                      FilterMask = filter2( Gauss(10,2),FilterMask); % proper mask
                        %% Filtration 2D
                        [Filtered_FFTData] = propFFT2_Filtering(FFTData,FilterMask);
                        clear FFTData NFrames thesh %FilterMask

                        %% IFFT2
                        [IFFT2DATA] = propIFFT2(Filtered_FFTData);
                        clear Filtered_FFTData
                        clear EData
                        
                        RMSF = abs(sqrt(sum(IFFT2DATA(:,:,1:512).^2,3)));
                        figure;
                        %surf(RMSF(2:end-1,2:end-1));shading interp; view(2); axis square;colorbar;
                        if(a>0)
                            surf(RMSF(a/2:nx+a/2-1,a/2:ny+a/2));shading interp; view(2); axis square;colorbar;
                        else
                           surf(RMSF);shading interp; view(2); axis square;colorbar;
                        end
                        caxis([0 0.1*max(max(RMSF))]); xlim([1 500]);ylim([1 500]);
                        colormap jet;
                        
                        %A=rms2image(RMS, filename);

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



