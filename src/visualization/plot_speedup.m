clear all; close all;
data_folder_path=['..',filesep,'..',filesep,'data',filesep,'raw',filesep,'num',filesep,'flat_shell',filesep];
folder = 'flat_shell_Jochen_wavefield_speedup_gpu_out';

% GPU
NDOFS = zeros(6,1);
Runtime_gpu = zeros(6,1);
nft = 4096*16;      % total number of samples
for k=1:6
    folder_name=[folder,filesep,num2str(k),'_output',filesep];
     filename_time=['time',num2str(k)'];
     filepath2=[data_folder_path,folder_name,filename_time];
    load(filepath2);% time vector t, TotalNofNodes, averageTime
    NDOFS(k,1) = TotalNofNodes*5/10^6;
    Runtime_gpu(k,1)= averageTime*nft/60;
end
figure;
plot(NDOFS,Runtime_gpu,'-o','Linewidth',2);
xlabel('NDOF *10^6');
ylabel('Run time [min]');
xlim([0.5 7]);
Runtime_gpu/60

% CPU
folder = 'flat_shell_Jochen_wavefield_speedup_cpu_out';
Runtime_cpu = zeros(6,1);
nft = 4096*16;      % total number of samples
for k=1:6
    folder_name=[folder,filesep,num2str(k),'_output',filesep];
     filename_time=['time',num2str(k)'];
     filepath2=[data_folder_path,folder_name,filename_time];
    load(filepath2);% time vector t, TotalNofNodes, averageTime
    NDOFS(k,1) = TotalNofNodes*5/10^6;
    Runtime_cpu(k,1)= averageTime*nft/60;
end
figure
plot(NDOFS,Runtime_cpu,'k-o','Linewidth',2);
xlabel('NDOF *10^6');
ylabel('Run time [min]');
xlim([0.5 7]);
Runtime_cpu/60

% speedup
speedup = Runtime_cpu./Runtime_gpu
figure;
fig_width = 6; fig_height = 5; 

plot(NDOFS,speedup,'k-o','Linewidth',1);
set(gca,'FontName','Times');
set(gcf,'Color','w');
set(gca,'Fontsize',10);
xlabel('NDOF *10^6','Fontsize',12);
ylabel('Speedup','Fontsize',12);
return;
%set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
% remove unnecessary white space
%set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
%fig.PaperPositionMode   = 'auto';
set(gcf,'PaperPositionMode','auto');
paper_fig_folder = 'E:\work\projects\nawa-bekker\ma-shm\reports\journal_papers\Elsevier\figs\';
figfilename='speedup';
print([paper_fig_folder,figfilename],'-dpng', '-r600'); 