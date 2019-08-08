clear all; close all;
% load wavefield data for the case of random excitation
% at actuator #1
% sensors amd defects placement according to Open Guided Waves platform (SHM plate)
% J. Moll , J. Kathol, C.-P. Fritzen, M. Moix-Bonet, M. Rennoch , M. Koerdt, A.S. Herrmann, M. GR Sause and M. Bach
% Open GuidedWaves: online platform for ultrasonic guided wave
% measurements, Structural Health Monitoring 1–12, 2018.
% see also: SHM datasets: Short instruction for working with the SHM data set

actuator_no=1;
fig_width=5;
fig_height=5;
data_folder_path=['..',filesep,'..',filesep,'data',filesep,'interim',filesep,'num',filesep,'flat_shell',filesep];
folder = 'flat_shell_Jochen_random_exc_signals';
caxis_cut=0.2;
% damage #5
for k=1:5
    folder_name=[folder,num2str(k),'_out',filesep,num2str(actuator_no),'_output',filesep];
    filename=['flat_shell_Vz_',num2str(actuator_no),'_500x500top'];
    filepath=[data_folder_path,folder_name,filename];
    load(filepath);% Data
    [m,n,nft]=size(Data);
    % RMS
    RMS=zeros(m,n);
    RMS_combined=zeros(m,n);
    for j=1:nft
        RMS=RMS+Data(:,:,j).^2;
    end
    Smax=max(max(RMS));
    figure;
    surf(RMS);shading interp;view(2);axis equal;colorbar;axis off;
    caxis([0 caxis_cut*Smax]);
    set(gcf,'Renderer','zbuffer');
    RMS_combined = RMS_combined + RMS;
end
Smax=max(max(RMS_combined(2:end-1,2:end-1)));
figure;
surf(RMS_combined(2:end-1,2:end-1));shading interp;view(2);axis equal;axis off;
caxis([0 caxis_cut*Smax]);
set(gcf,'Renderer','zbuffer');
set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
% remove unnecessary white space
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
set(gcf,'PaperPositionMode','auto');    
print('RMS5','-dpng', '-r600');
% damage #14
for k=6:10
    folder_name=[folder,num2str(k),'_out',filesep,num2str(actuator_no),'_output',filesep];
    filename=['flat_shell_Vz_',num2str(actuator_no),'_500x500top'];
    filepath=[data_folder_path,folder_name,filename];
    load(filepath);% Data
    [m,n,nft]=size(Data);
    % RMS
    RMS=zeros(m,n);
    RMS_combined=zeros(m,n);
    for j=1:nft
        RMS=RMS+Data(:,:,j).^2;
    end
    Smax=max(max(RMS));
    figure;
    surf(RMS);shading interp;view(2);axis equal;colorbar;axis off;
    caxis([0 caxis_cut*Smax]);
    set(gcf,'Renderer','zbuffer');
    RMS_combined = RMS_combined + RMS;
end
Smax=max(max(RMS(2:end-1,2:end-1)));
figure;
surf(RMS(2:end-1,2:end-1));shading interp;view(2);axis equal;axis off;
caxis([0 caxis_cut*Smax]);
set(gcf,'Renderer','zbuffer');
set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
% remove unnecessary white space
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
set(gcf,'PaperPositionMode','auto');    
print('RMS14','-dpng', '-r600');
% damage #22
for k=11:15
    folder_name=[folder,num2str(k),'_out',filesep,num2str(actuator_no),'_output',filesep];
    filename=['flat_shell_Vz_',num2str(actuator_no),'_500x500top'];
    filepath=[data_folder_path,folder_name,filename];
    load(filepath);% Data
    [m,n,nft]=size(Data);
    % RMS
    RMS=zeros(m,n);
    RMS_combined=zeros(m,n);
    for j=1:nft
        RMS=RMS+Data(:,:,j).^2;
    end
    Smax=max(max(RMS));
    figure;
    surf(RMS);shading interp;view(2);axis equal;colorbar;axis off;
    caxis([0 caxis_cut*Smax]);
    set(gcf,'Renderer','zbuffer');
    RMS_combined = RMS_combined + RMS;
end
Smax=max(max(RMS(2:end-1,2:end-1)));
figure;
surf(RMS(2:end-1,2:end-1));shading interp;view(2);axis equal;axis off;
caxis([0 caxis_cut*Smax]);
set(gcf,'Renderer','zbuffer');
set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
% remove unnecessary white space
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
set(gcf,'PaperPositionMode','auto');    
print('RMS22','-dpng', '-r600');
%% wavenumber filtering
Nmed = 1; k = 1; thresh_scaling = 1; NFrames = 1:128;%NFrames = 1:64;
PRT = 1;% if PRT = 1 then print
blur = 1;   % guassian blur of the filter mask %[ 0.5, 1 lub 2]
ct = 0;     %cut
%[H2D] = Hann2D(size(Data,2),size(Data,1),0.1); %st 0.1 %tramp 0.00001
% H2D = ones(size(Data,1),size(Data,2));
% 
% %% Color map
% whitethreshold = .05;
% blackthreshold = .05;
% [CmapB] = colormapwlasny(whitethreshold,blackthreshold);
% Cmap = jet(256);  

%% median filtering
% [FData] = propMedian(Data,Nmed);
%clear Data

%% FFT 2D
[FFTData] = propFFT2(Data);
%% filter mask
[FilterMask] = propMask(FFTData,NFrames,thresh_scaling);
%% Filtation 2D
[Filtered_FFTData] = propFFT2_Filtering(FFTData,FilterMask);
clear FFTData NFrames thesh_scaling FilterMask
%% IFFT2
[IFFT2DATA] = propIFFT2(Filtered_FFTData);
clear Filtered_FFTData
%% RMSF
caxis_cut=0.5;
[RMSF] = propRMS(IFFT2DATA);
Smax=max(max(RMSF(2:end-1,2:end-1)));
figure;
surf(RMSF(2:end-1,2:end-1));shading interp;view(2);axis equal;axis off;
caxis([0 caxis_cut*Smax]);
set(gcf,'Renderer','zbuffer');
set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
% remove unnecessary white space
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
set(gcf,'PaperPositionMode','auto');    
print('RMSF22','-dpng', '-r600');