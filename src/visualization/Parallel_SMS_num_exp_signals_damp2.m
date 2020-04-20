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
% prepare figure output path
modelname = name; 
figure_output_path = prepare_figure_paths(modelname);
modelfolder = 'flat_shell'; % name of folder
modelname =  'flat_shell_Parallel_SMS';
% prepare model input/output paths
model_interim_path = prepare_model_paths('interim','num',modelfolder,modelname);

%% Input for flat_shell
input_file_no = [21,22,23,24,25,26,27,28,29,30,31,32]; % frequency = [16.5e3/5,50e3/5,100e3/5];
%% input for experiment
list = {'491x491p_16_5kHz_5HC_x5_15Vpp', ...     
         '492x492p_50kHz_5HC_x20_15Vpp', ...
         '491x491p_100kHz_5HC_x20_15Vpp'};  
% create path to the experimental data folder
input_data_path = fullfile( projectroot, 'data','raw','exp', filesep );
start_index = 10; 
%% input for figures
fig_width =12; % figure widht in cm
fig_height=10; % figure height in cm
lc='rgbm';
% selected points for comparison
% coordinates C(0.3,0.4506), B(0.3,0.3), A(0.4506,0.3), D(0.1506,0.1506)
ind_Ax=245; ind_Ay=368;
ind_Bx=245; ind_By=245;
ind_Cx=368; ind_Cy=245;
ind_Dx=123; ind_Dy=123;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 16.5 kHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical
tasks=[4];% 16.5 kHz
% experiment
filename = list{1}; %16.5 kHz
load([input_data_path,filename]); % Data, time
[nx,ny,nt]=size(Data);
Data_exp=Data(1:nx,1:ny,:);
sA_exp = squeeze(Data_exp(ind_Ax,ind_Ay,:));
sB_exp = squeeze(Data_exp(ind_Bx,ind_By,:));
sC_exp = squeeze(Data_exp(ind_Cx,ind_Cy,:));
sD_exp = squeeze(Data_exp(ind_Dx,ind_Dy,:));
% numerical damped
sA_num=zeros(nt,length(tasks)+1)+NaN;
sB_num=zeros(nt,length(tasks)+1)+NaN;
sC_num=zeros(nt,length(tasks)+1)+NaN;
sD_num=zeros(nt,length(tasks)+1)+NaN;
t=zeros(nt,length(tasks)+1);
alphaz=zeros(length(tasks),1);
alphaxy=zeros(length(tasks),1);
%deltat= 7.8*1e-6;
deltat= 0;
c=0;
for test_case = tasks
    c=c+1;
    input_no = input_file_no(test_case);
    run(fullfile('inputs',['input',num2str(input_no)]));
    alphaz(c) = etad_z;
    alphaxy(c) = etad_xy;
    interim_output_name = [model_interim_path,filesep,num2str(test_case),'_output',filesep];
    if(overwrite)
        load([interim_output_name,'flat_shell_Vz_',num2str(test_case),'_982x982bottom']);
        Data_num=Data(1:nx,1:ny,:);% Data
        [~,~,nt_num]=size(Data_num);
        load([interim_output_name,'t_frames']); % t_frames (time vector)
        t(1:length(t_frames),c) = t_frames;  
        % remove deltat delay
        [Data_num] = wavefield_delta_t_correction(-1*Data_num,t_frames,deltat,1);
        % calculate amplification factor
        if(c==1)
            [Am,A] = amplification_factor(Data_exp(start_index:nx-start_index*4,start_index:ny-start_index*4,1:nt_num),Data_num(start_index:nx-start_index*4,start_index:ny-start_index*4,1:nt_num));
        end
        sA_num(1:nt_num,c) = Am*Data_num(ind_Ax,ind_Ay,:);
        sB_num(1:nt_num,c) = Am*Data_num(ind_Bx,ind_By,:);
        sC_num(1:nt_num,c) = Am*Data_num(ind_Cx,ind_Cy,:);
        sD_num(1:nt_num,c) = Am*Data_num(ind_Dx,ind_Dy,:);
    else
        fprintf([modelname,' test case: %d already exist\n'], test_case);
    end
end
sA_num = (-1)*sA_num*1e3; % scale to mm/s
sB_num = sB_num*1e3; % scale to mm/s
sC_num = sC_num*1e3; % scale to mm/s
sD_num = sD_num*1e3; % scale to mm/s
sA_exp = sA_exp*1e3; % scale to mm/s
sB_exp = sB_exp*1e3; % scale to mm/s
sC_exp = sC_exp*1e3; % scale to mm/s
sD_exp = sD_exp*1e3; % scale to mm/s
time = time*1e6; % scale to micro seconds
t = t*1e6; % scale to micro seconds
%% POINT A AND B 16.5 kHz
hp=cell(length(tasks)+1,1);
figure;
subplot(2,2,1);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sA_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sA_exp,'k:','LineWidth',2);
Smax=1.1*max(max(abs(sA_exp)));
Smin=-Smax;
axis([0 time(256) Smin Smax]);
set(gcf,'color','white');
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
ylabel('Amplitude [mm/s]','Fontsize',11);
title(['A, ', num2str(f_2/1000),' [kHz]'],'Fontsize',11); 
subplot(2,2,2);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sB_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sB_exp,'k:','LineWidth',2);
axis([0 time(256) Smin Smax]);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
ylabel('Amplitude [mm/s]','Fontsize',11);
title(['B, ', num2str(f_2/1000),' [kHz]'],'Fontsize',11); 
subplot(2,2,3);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sA_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sA_exp,'k:','LineWidth',2);
%lg=legend([hp{1},hp{2},hp{3},hp{4}],{['$\alpha =$ ',num2str(alphaz(1))],['$\alpha =$ ',num2str(alphaz(2))],['$\alpha =$ ',num2str(alphaz(3))],'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
lg=legend([hp{1}],{['$\alpha =$ ',num2str(alphaz(1)),' num.']},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
legend boxoff;
xlabel({['t [','$\mu$','s]']},'Fontsize',11,'interpreter','latex');
ylabel('Amplitude [mm/s]','Fontsize',11);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
axis([400 500 Smin Smax]);
subplot(2,2,4);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sB_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sB_exp,'k:','LineWidth',2);
%hLeg=legend([hp{1},hp{2},hp{3},hp{4}],{['$\alpha =$ ',num2str(alphaz(1))],['$\alpha =$ ',num2str(alphaz(2))],['$\alpha =$ ',num2str(alphaz(3))],'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
hLeg=legend([hp{2}],{'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
%set(hLeg,'visible','off')
legend boxoff;
axis([500 600 Smin Smax]);
xlabel({['t [','$\mu$','s]']},'Fontsize',11,'interpreter','latex');
ylabel('Amplitude [mm/s]','Fontsize',11);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
set(gcf, 'Units','centimeters', 'Position',[4 4 fig_width fig_height]);
figfilename = ['num_exp_signals_damp_A_B_',num2str(f_2/1000),'_kHz.png'];
print([figure_output_path,figfilename],'-dpng', '-r600'); 
%% POINT C AND D 16.5 kHz
figure;
subplot(2,2,1);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sC_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sC_exp,'k:','LineWidth',2);
Smax=1.1*max(max(abs(sC_exp)));
Smin=-Smax;
axis([0 time(256) Smin Smax]);
set(gcf,'color','white');
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
ylabel('Amplitude [mm/s]','Fontsize',11);
title(['C, ', num2str(f_2/1000),' [kHz]'],'Fontsize',11); 
subplot(2,2,2);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sD_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sD_exp,'k:','LineWidth',2);
axis([0 time(256) Smin Smax]);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
ylabel('Amplitude [mm/s]','Fontsize',11);
title(['D, ', num2str(f_2/1000),' [kHz]'],'Fontsize',11); 
subplot(2,2,3);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sC_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sC_exp,'k:','LineWidth',2);
%lg=legend([hp{1},hp{2},hp{3},hp{4}],{['$\alpha =$ ',num2str(alphaz(1))],['$\alpha =$ ',num2str(alphaz(2))],['$\alpha =$ ',num2str(alphaz(3))],'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
lg=legend([hp{1}],{['$\alpha =$ ',num2str(alphaz(1)),' num.']},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
legend boxoff;
xlabel({['t [','$\mu$','s]']},'Fontsize',11,'interpreter','latex');
ylabel('Amplitude [mm/s]','Fontsize',11);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
axis([400 500 Smin Smax]);
subplot(2,2,4);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sD_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sD_exp,'k:','LineWidth',2);
%hLeg=legend([hp{1},hp{2},hp{3},hp{4}],{['$\alpha =$ ',num2str(alphaz(1))],['$\alpha =$ ',num2str(alphaz(2))],['$\alpha =$ ',num2str(alphaz(3))],'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
hLeg=legend([hp{2}],{'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
%set(hLeg,'visible','off')
legend boxoff;
axis([700 800 Smin Smax]);
xlabel({['t [','$\mu$','s]']},'Fontsize',11,'interpreter','latex');
ylabel('Amplitude [mm/s]','Fontsize',11);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
set(gcf, 'Units','centimeters', 'Position',[4 4 fig_width fig_height]);
figfilename = ['num_exp_signals_damp_C_D_',num2str(f_2/1000),'_kHz.png'];
print([figure_output_path,figfilename],'-dpng', '-r600'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 50 kHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical
tasks=[5];% 50 kHz
% experiment
filename = list{2}; %50 kHz
load([input_data_path,filename]); % Data, time
[nx,ny,nt]=size(Data);
Data_exp=Data(1:nx,1:ny,:);
sA_exp = squeeze(Data_exp(ind_Ax,ind_Ay,:));
sB_exp = squeeze(Data_exp(ind_Bx,ind_By,:));
sC_exp = squeeze(Data_exp(ind_Cx,ind_Cy,:));
sD_exp = squeeze(Data_exp(ind_Dx,ind_Dy,:));
% numerical damped
sA_num=zeros(nt,length(tasks)+1)+NaN;
sB_num=zeros(nt,length(tasks)+1)+NaN;
sC_num=zeros(nt,length(tasks)+1)+NaN;
sD_num=zeros(nt,length(tasks)+1)+NaN;
t=zeros(nt,length(tasks)+1);
alphaz=zeros(length(tasks),1);
alphaxy=zeros(length(tasks),1);
%deltat= 7.8*1e-6;
deltat= 0;
c=0;
for test_case = tasks
    c=c+1;
    input_no = input_file_no(test_case);
    run(fullfile('inputs',['input',num2str(input_no)]));
    alphaz(c) = etad_z;
    alphaxy(c) = etad_xy;
    interim_output_name = [model_interim_path,filesep,num2str(test_case),'_output',filesep];
    if(overwrite)
        load([interim_output_name,'flat_shell_Vz_',num2str(test_case),'_982x982bottom']);
        Data_num=Data(1:nx,1:ny,:);% Data
        [~,~,nt_num]=size(Data_num);
        load([interim_output_name,'t_frames']); % t_frames (time vector)
        t(1:length(t_frames),c) = t_frames;  
        % remove deltat delay
        [Data_num] = wavefield_delta_t_correction(-1*Data_num,t_frames,deltat,1);
        % calculate amplification factor
        if(c==1)
            [Am,A] = amplification_factor(Data_exp(start_index:nx-start_index*4,start_index:ny-start_index*4,1:nt_num),Data_num(start_index:nx-start_index*4,start_index:ny-start_index*4,1:nt_num));
        end
        sA_num(1:nt_num,c) = Am*Data_num(ind_Ax,ind_Ay,:);
        sB_num(1:nt_num,c) = Am*Data_num(ind_Bx,ind_By,:);
        sC_num(1:nt_num,c) = Am*Data_num(ind_Cx,ind_Cy,:);
        sD_num(1:nt_num,c) = Am*Data_num(ind_Dx,ind_Dy,:);
    else
        fprintf([modelname,' test case: %d already exist\n'], test_case);
    end
end
sA_num = (-1)*sA_num*1e3; % scale to mm/s
sB_num = sB_num*1e3; % scale to mm/s
sC_num = sC_num*1e3; % scale to mm/s
sD_num = (-1)*sD_num*1e3; % scale to mm/s
sA_exp = sA_exp*1e3; % scale to mm/s
sB_exp = sB_exp*1e3; % scale to mm/s
sC_exp = sC_exp*1e3; % scale to mm/s
sD_exp = sD_exp*1e3; % scale to mm/s
time = time*1e6; % scale to micro seconds
t = t*1e6; % scale to micro seconds
%% POINT A AND B 50 kHz
hp=cell(length(tasks)+1,1);
figure;
subplot(2,2,1);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sA_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sA_exp,'k:','LineWidth',2);
Smax=1.1*max(max(abs(sA_exp)));
Smin=-Smax;
axis([0 time(256) Smin Smax]);
set(gcf,'color','white');
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
ylabel('Amplitude [mm/s]','Fontsize',11);
title(['A, ', num2str(f_2/1000),' [kHz]'],'Fontsize',11); 
subplot(2,2,2);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sB_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sB_exp,'k:','LineWidth',2);
axis([0 time(256) Smin Smax]);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
ylabel('Amplitude [mm/s]','Fontsize',11);
title(['B, ', num2str(f_2/1000),' [kHz]'],'Fontsize',11); 
subplot(2,2,3);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sA_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sA_exp,'k:','LineWidth',2);
%lg=legend([hp{1},hp{2},hp{3},hp{4}],{['$\alpha =$ ',num2str(alphaz(1))],['$\alpha =$ ',num2str(alphaz(2))],['$\alpha =$ ',num2str(alphaz(3))],'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
lg=legend([hp{1}],{['$\alpha =$ ',num2str(alphaz(1)),' num.']},'Fontsize',9,'interpreter','latex','Location','NorthOutside');

legend boxoff;
xlabel({['t [','$\mu$','s]']},'Fontsize',11,'interpreter','latex');
ylabel('Amplitude [mm/s]','Fontsize',11);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
axis([250 280 Smin Smax]);
subplot(2,2,4);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sB_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sB_exp,'k:','LineWidth',2);
%hLeg=legend([hp{1},hp{2},hp{3},hp{4}],{['$\alpha =$ ',num2str(alphaz(1))],['$\alpha =$ ',num2str(alphaz(2))],['$\alpha =$ ',num2str(alphaz(3))],'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
hLeg=legend([hp{2}],{'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
%set(hLeg,'visible','off')
legend boxoff;
axis([350 380 Smin Smax]);
xlabel({['t [','$\mu$','s]']},'Fontsize',11,'interpreter','latex');
ylabel('Amplitude [mm/s]','Fontsize',11);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
set(gcf, 'Units','centimeters', 'Position',[4 4 fig_width fig_height]);
figfilename = ['num_exp_signals_damp_A_B_',num2str(f_2/1000),'_kHz.png'];
print([figure_output_path,figfilename],'-dpng', '-r600'); 
%% POINT C AND D 50 kHz
figure;
subplot(2,2,1);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sC_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sC_exp,'k:','LineWidth',2);
Smax=1.1*max(max(abs(sC_exp)));
Smin=-Smax;
axis([0 time(384) Smin Smax]);
set(gcf,'color','white');
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
ylabel('Amplitude [mm/s]','Fontsize',11);
title(['C, ', num2str(f_2/1000),' [kHz]'],'Fontsize',11); 
subplot(2,2,2);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sD_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sD_exp,'k:','LineWidth',2);
axis([0 time(384) 0.5*Smin 0.5*Smax]);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
ylabel('Amplitude [mm/s]','Fontsize',11);
title(['D, ', num2str(f_2/1000),' [kHz]'],'Fontsize',11); 
subplot(2,2,3);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sC_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sC_exp,'k:','LineWidth',2);
%lg=legend([hp{1},hp{2},hp{3},hp{4}],{['$\alpha =$ ',num2str(alphaz(1))],['$\alpha =$ ',num2str(alphaz(2))],['$\alpha =$ ',num2str(alphaz(3))],'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
lg=legend([hp{1}],{['$\alpha =$ ',num2str(alphaz(1)),' num.']},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
legend boxoff;
xlabel({['t [','$\mu$','s]']},'Fontsize',11,'interpreter','latex');
ylabel('Amplitude [mm/s]','Fontsize',11);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
axis([300 320 Smin Smax]);
subplot(2,2,4);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sD_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sD_exp,'k:','LineWidth',2);
%hLeg=legend([hp{1},hp{2},hp{3},hp{4}],{['$\alpha =$ ',num2str(alphaz(1))],['$\alpha =$ ',num2str(alphaz(2))],['$\alpha =$ ',num2str(alphaz(3))],'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
hLeg=legend([hp{2}],{'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
%set(hLeg,'visible','off')
legend boxoff;
axis([500 520 0.5*Smin 0.5*Smax]);
xlabel({['t [','$\mu$','s]']},'Fontsize',11,'interpreter','latex');
ylabel('Amplitude [mm/s]','Fontsize',11);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
set(gcf, 'Units','centimeters', 'Position',[4 4 fig_width fig_height]);
figfilename = ['num_exp_signals_damp_C_D_',num2str(f_2/1000),'_kHz.png'];
print([figure_output_path,figfilename],'-dpng', '-r600'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 100 kHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical
tasks=[6];% 100 kHz
% experiment
filename = list{3}; %100 kHz
load([input_data_path,filename]); % Data, time
[nx,ny,nt]=size(Data);
Data_exp=Data(1:nx,1:ny,:);
sA_exp = squeeze(Data_exp(ind_Ax,ind_Ay,:));
sB_exp = squeeze(Data_exp(ind_Bx,ind_By,:));
sC_exp = squeeze(Data_exp(ind_Cx,ind_Cy,:));
sD_exp = squeeze(Data_exp(ind_Dx,ind_Dy,:));
% numerical damped
sA_num=zeros(nt,length(tasks)+1)+NaN;
sB_num=zeros(nt,length(tasks)+1)+NaN;
sC_num=zeros(nt,length(tasks)+1)+NaN;
sD_num=zeros(nt,length(tasks)+1)+NaN;
t=zeros(nt,length(tasks)+1);
alphaz=zeros(length(tasks),1);
alphaxy=zeros(length(tasks),1);
%deltat= 7.8*1e-6;
deltat= 0;
c=0;
for test_case = tasks
    c=c+1;
    input_no = input_file_no(test_case);
    run(fullfile('inputs',['input',num2str(input_no)]));
    alphaz(c) = etad_z;
    alphaxy(c) = etad_xy;
    interim_output_name = [model_interim_path,filesep,num2str(test_case),'_output',filesep];
    if(overwrite)
        load([interim_output_name,'flat_shell_Vz_',num2str(test_case),'_982x982bottom']);
        Data_num=Data(1:nx,1:ny,:);% Data
        [~,~,nt_num]=size(Data_num);
        load([interim_output_name,'t_frames']); % t_frames (time vector)
        t(1:length(t_frames),c) = t_frames;  
        % remove deltat delay
        [Data_num] = wavefield_delta_t_correction(-1*Data_num,t_frames,deltat,1);
        % calculate amplification factor
        if(c==1)
            [Am,A] = amplification_factor(Data_exp(start_index:nx-start_index*4,start_index:ny-start_index*4,1:nt_num),Data_num(start_index:nx-start_index*4,start_index:ny-start_index*4,1:nt_num));
        end
        sA_num(1:nt_num,c) = Am*Data_num(ind_Ax,ind_Ay,:);
        sB_num(1:nt_num,c) = Am*Data_num(ind_Bx,ind_By,:);
        sC_num(1:nt_num,c) = Am*Data_num(ind_Cx,ind_Cy,:);
        sD_num(1:nt_num,c) = Am*Data_num(ind_Dx,ind_Dy,:);
    else
        fprintf([modelname,' test case: %d already exist\n'], test_case);
    end
end
sA_num = (-1)*sA_num*1e3; % scale to mm/s
sB_num = sB_num*1e3; % scale to mm/s
sC_num = sC_num*1e3; % scale to mm/s
sD_num = (-1)*sD_num*1e3; % scale to mm/s
sA_exp = sA_exp*1e3; % scale to mm/s
sB_exp = sB_exp*1e3; % scale to mm/s
sC_exp = sC_exp*1e3; % scale to mm/s
sD_exp = sD_exp*1e3; % scale to mm/s
time = time*1e6; % scale to micro seconds
t = t*1e6; % scale to micro seconds
%% POINT A AND B 100 kHz
hp=cell(length(tasks)+1,1);
figure;
subplot(2,2,1);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sA_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sA_exp,'k:','LineWidth',2);
Smax=1.1*max(max(abs(sA_exp)));
Smin=-Smax;
axis([0 350 Smin Smax]);
set(gcf,'color','white');
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
ylabel('Amplitude [mm/s]','Fontsize',11);
title(['A, ', num2str(f_2/1000),' [kHz]'],'Fontsize',11); 
subplot(2,2,2);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sB_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sB_exp,'k:','LineWidth',2);
axis([0 350 Smin Smax]);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
ylabel('Amplitude [mm/s]','Fontsize',11);
title(['B, ', num2str(f_2/1000),' [kHz]'],'Fontsize',11); 
subplot(2,2,3);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sA_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sA_exp,'k:','LineWidth',2);
%lg=legend([hp{1},hp{2},hp{3},hp{4}],{['$\alpha =$ ',num2str(alphaz(1))],['$\alpha =$ ',num2str(alphaz(2))],['$\alpha =$ ',num2str(alphaz(3))],'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
lg=legend([hp{1}],{['$\alpha =$ ',num2str(alphaz(1)),' num.']},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
legend boxoff;
xlabel({['t [','$\mu$','s]']},'Fontsize',11,'interpreter','latex');
ylabel('Amplitude [mm/s]','Fontsize',11);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
axis([220 235 Smin Smax]);
subplot(2,2,4);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sB_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sB_exp,'k:','LineWidth',2);
%hLeg=legend([hp{1},hp{2},hp{3},hp{4}],{['$\alpha =$ ',num2str(alphaz(1))],['$\alpha =$ ',num2str(alphaz(2))],['$\alpha =$ ',num2str(alphaz(3))],'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
hLeg=legend([hp{2}],{'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
%set(hLeg,'visible','off')
legend boxoff;
axis([300 315 Smin Smax]);
xlabel({['t [','$\mu$','s]']},'Fontsize',11,'interpreter','latex');
ylabel('Amplitude [mm/s]','Fontsize',11);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
set(gcf, 'Units','centimeters', 'Position',[4 4 fig_width fig_height]);
figfilename = ['num_exp_signals_damp_A_B_',num2str(f_2/1000),'_kHz.png'];
print([figure_output_path,figfilename],'-dpng', '-r600'); 
%% POINT C AND D 100 kHz
figure;
subplot(2,2,1);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sC_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sC_exp,'k:','LineWidth',2);
Smax=1.1*max(max(abs(sC_exp)));
Smin=-Smax;
axis([0 480 Smin Smax]);
set(gcf,'color','white');
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
ylabel('Amplitude [mm/s]','Fontsize',11);
title(['C, ', num2str(f_2/1000),' [kHz]'],'Fontsize',11); 
subplot(2,2,2);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sD_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sD_exp,'k:','LineWidth',2);
axis([0 480 0.2*Smin 0.2*Smax]);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
ylabel('Amplitude [mm/s]','Fontsize',11);
title(['D, ', num2str(f_2/1000),' [kHz]'],'Fontsize',11); 
subplot(2,2,3);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sC_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sC_exp,'k:','LineWidth',2);
%lg=legend([hp{1},hp{2},hp{3},hp{4}],{['$\alpha =$ ',num2str(alphaz(1))],['$\alpha =$ ',num2str(alphaz(2))],['$\alpha =$ ',num2str(alphaz(3))],'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
lg=legend([hp{1}],{['$\alpha =$ ',num2str(alphaz(1)),' num.']},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
legend boxoff;
xlabel({['t [','$\mu$','s]']},'Fontsize',11,'interpreter','latex');
ylabel('Amplitude [mm/s]','Fontsize',11);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
axis([250 265 Smin Smax]);
subplot(2,2,4);
for c=1:length(tasks)
    hp{c}=plot((t(:,c)-deltat),sD_num(:,c) ,lc(c),'LineWidth',1);
    hold on;
end
hp{c+1}=plot(time,sD_exp,'k:','LineWidth',2);
%hLeg=legend([hp{1},hp{2},hp{3},hp{4}],{['$\alpha =$ ',num2str(alphaz(1))],['$\alpha =$ ',num2str(alphaz(2))],['$\alpha =$ ',num2str(alphaz(3))],'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
hLeg=legend([hp{2}],{'experiment'},'Fontsize',9,'interpreter','latex','Location','NorthOutside');
%set(hLeg,'visible','off')
legend boxoff;
axis([440 455 0.2*Smin 0.2*Smax]);
xlabel({['t [','$\mu$','s]']},'Fontsize',11,'interpreter','latex');
ylabel('Amplitude [mm/s]','Fontsize',11);
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
set(gcf, 'Units','centimeters', 'Position',[4 4 fig_width fig_height]);
figfilename = ['num_exp_signals_damp_C_D_',num2str(f_2/1000),'_kHz.png'];
print([figure_output_path,figfilename],'-dpng', '-r600'); 

