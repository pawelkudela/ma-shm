clear all; close all;
% load wavefield data for the case of random excitation
% at actuator #1
% sensors amd defects placement according to Open Guided Waves dataset (SHM plate)
% J. Moll , J. Kathol, C.-P. Fritzen, M. Moix-Bonet, M. Rennoch , M. Koerdt, A.S. Herrmann, M. GR Sause and M. Bach
% Open GuidedWaves: online platform for ultrasonic guided wave
% measurements, Structural Health Monitoring 1–12, 2018.
% see also: SHM datasets: Short instruction for working with the SHM data set
% E:\work\projects\nawa-bekker\ma-shm\data\raw\num\flat_shell\flat_shell_Jochen_random_exc_signals1_out\1_output
% ..\..\data\raw\num\flat_shell\flat_shell_Jochen_random_exc_signals1_out\1_output\voltage_1':
actuator_no=1;
data_folder_path=['..',filesep,'..',filesep,'data',filesep,'raw',filesep,'num',filesep,'flat_shell',filesep];
folder = 'flat_shell_Jochen_random_exc_signals';
nft=122880;
%% damage #5
% sensor signal data (number of time steps, sensor no, random excitation number)
Data5=zeros(nft,12,5); 
% random excitation signal
exc_sig5=zeros(nft,5); 
% time vector t is the same across all cases
for k=1:5
    folder_name=[folder,num2str(k),'_out',filesep,num2str(actuator_no),'_output',filesep];
    filename_voltage=['voltage',num2str(actuator_no)','.mat'];
    filename_time=['time',num2str(actuator_no)'];
    filepath1=[data_folder_path,folder_name,filename_voltage];
    load(filepath1);% voltage
    filepath2=[data_folder_path,folder_name,filename_time];
    load(filepath2);% time vector t and excitation signal st
    Data5(:,:,k) = voltage;
    exc_sig5(:,k) = st; 
end
% plot selected data
% plot random excitation signals
fig_width = 12; fig_height = 5; 
linewidth = 0.5;
figure;
set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height],'Color','w'); 
hold on;
plot(t*1e6,exc_sig5,'LineWidth',2*linewidth);
set(gca,'FontName','Times');
legend('exc1','exc2','exc3','exc4','exc5');
xlabel('Time [µs]');
ylabel('Amplitude [V]');
xlim([0 1200]);
box on;
paper_fig_folder = 'E:\work\projects\nawa-bekker\ma-shm\reports\journal_papers\Elsevier\figs\';
figure;
set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height],'Color','w'); 
hold on;
plot(t*1e6,exc_sig5(:,5),'k-','LineWidth',2*linewidth);
set(gca,'FontName','Times');
title('Random excitation signal');
xlabel('Time [µs]');
ylabel('Amplitude [V]');
xlim([0 1200]);
box on;
figname='random_exc_sig';
figfilename=fullfile(paper_fig_folder,'beamer_figs',figname);
print(figfilename,'-dpng', '-r600'); 

% plot signal registerd in sensor 7 for random excitation 2
figure;
hold on;
plot(t,Data5(:,7,2));
% plot signal registerd in sensor 10 for random excitation 2
plot(t,Data5(:,10,2),'r');
legend('sensor 7','sensor 10');

%% damage #14
% sensor signal data (number of time steps, sensor no, random excitation number)
Data14=zeros(nft,12,5); 
% random excitation signal
exc_sig14=zeros(nft,5); 
% time vector t is the same across all cases
for k=6:10
    folder_name=[folder,num2str(k),'_out',filesep,num2str(actuator_no),'_output',filesep];
    filename_voltage=['voltage',num2str(actuator_no)','.mat'];
    filename_time=['time',num2str(actuator_no)'];
    filepath1=[data_folder_path,folder_name,filename_voltage];
    load(filepath1);% voltage
    filepath2=[data_folder_path,folder_name,filename_time];
    load(filepath2);% time vector t and excitation signal st
    Data14(:,:,k) = voltage;
    exc_sig14(:,k) = st; 
end

%% damage #22
% sensor signal data (number of time steps, sensor no, random excitation number)
Data22=zeros(nft,12,5); 
% random excitation signal
exc_sig22=zeros(nft,5); 
% time vector t is the same across all cases
for k=11:15
    folder_name=[folder,num2str(k),'_out',filesep,num2str(actuator_no),'_output',filesep];
    filename_voltage=['voltage',num2str(actuator_no)','.mat'];
    filename_time=['time',num2str(actuator_no)'];
    filepath1=[data_folder_path,folder_name,filename_voltage];
    load(filepath1);% voltage
    filepath2=[data_folder_path,folder_name,filename_time];
    load(filepath2);% time vector t and excitation signal st
    Data22(:,:,k) = voltage;
    exc_sig22(:,k) = st; 
end