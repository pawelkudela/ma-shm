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
modelname =  'flat_shell_EWSHM2020';
% prepare model input/output paths
model_interim_path = prepare_model_paths('interim','num',modelfolder,modelname);

%% Input for flat_shell
input_file_no = [21,22,23,24,25,26,27,28,29,30,31,32,33,34,35]; % frequency = [16.5e3/5,50e3/5,100e3/5];
%% input for experiment
list = {'491x491p_16_5kHz_5HC_x5_15Vpp', ...     
         '492x492p_50kHz_5HC_x20_15Vpp', ...
         '491x491p_100kHz_5HC_x20_15Vpp'};  
% create path to the experimental data folder
input_data_path = fullfile( projectroot, 'data','raw','exp', filesep );
%% input for figures
fig_width =5; % figure widht in cm
fig_height=5; % figure height in cm
lc='rgbm';

%% 16.5 kHz
% numerical
tasks=[4,7,10,13];% 16.5 kHz
c=0;
alphaz=zeros(length(tasks),1);
alphaxy=zeros(length(tasks),1);
E=zeros(1024,length(tasks)+1)+NaN;
t=zeros(1024,length(tasks)+1);
for test_case = tasks
    c=c+1;
    input_no = input_file_no(test_case);
    run(fullfile('inputs',['input',num2str(input_no)]));
    alphaz(c) = etad_z;
    alphaxy(c) = etad_xy;
    interim_output_name = [model_interim_path,filesep,num2str(test_case),'_output',filesep];
    if(overwrite)
            load([interim_output_name,'flat_shell_Vz_',num2str(test_case),'_982x982bottom']);
            Data=6*Data(1:491,1:491,:);% Data
            load([interim_output_name,'t_frames']); % t_frames (time vector)

           [En] = wavefield_energy(Data);
           E(1:length(En),c) = En;  
           t(1:length(t_frames),c) = t_frames;  
    else
          fprintf([modelname,' test case: %d already exist\n'], test_case);
    end
end

% experiment
filename = list{1}; %16.5 kHz
load([input_data_path,filename]); % Data, time
[En] = wavefield_energy(Data);
E(1:length(En),c+1) = En;  
deltat= 7.8*1e-6;
for c=1:length(tasks)
    plot((t(:,c)+deltat)*1e3,E(:,c),lc(c),'LineWidth',1);
    hold on;
end
plot(time*1e3,En,'k','LineWidth',2);
axis([0 t_frames(end)*1e3 0 1]);

set(gcf,'color','white');
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
xlabel('t [ms]','Fontsize',12);
ylabel('E [-]','Fontsize',12);
hold on;
title(['Frequency ',num2str(f_2/1000),' [kHz]'],'Fontsize',12); 
legend({['$\alpha_{z}=$ ',num2str(alphaz(1)),', $\alpha_{xy}=$ ',num2str(alphaxy(1))],['$\alpha_{z}=$ ',num2str(alphaz(2)),', $\alpha_{xy}=$ ',num2str(alphaxy(2))],['$\alpha_{z}=$ ',num2str(alphaz(3)),', $\alpha_{xy}=$ ',num2str(alphaxy(3))],['$\alpha_{z}=$ ',num2str(alphaz(4)),', $\alpha_{xy}=$ ',num2str(alphaxy(4))],'experiment'},'Fontsize',12,'interpreter','latex');
figfilename = ['num_exp_energy',num2str(f_2/1000),'_kHz.png'];
print([figure_output_path,figfilename],'-dpng', '-r600'); 
close all;
%% 50 kHz
% numerical
tasks=[5,8,11,14];% 50 kHz
c=0;
alphaz=zeros(length(tasks),1);
alphaxy=zeros(length(tasks),1);
E=zeros(1024,length(tasks)+1)+NaN;
t=zeros(1024,length(tasks)+1);
for test_case = tasks
    c=c+1;
    input_no = input_file_no(test_case);
    run(fullfile('inputs',['input',num2str(input_no)]));
    alphaz(c) = etad_z;
    alphaxy(c) = etad_xy;
    interim_output_name = [model_interim_path,filesep,num2str(test_case),'_output',filesep];

      if(overwrite)
            load([interim_output_name,'flat_shell_Vz_',num2str(test_case),'_982x982bottom']);
            Data=6*Data(1:491,1:491,:);% Data
            load([interim_output_name,'t_frames']); % t_frames (time vector)

           [En] = wavefield_energy(Data);
           E(1:length(En),c) = En;   
           t(1:length(t_frames),c) = t_frames;
      else
          fprintf([modelname,' test case: %d already exist\n'], test_case);
      end
end

% experiment
filename = list{2}; %50 kHz
load([input_data_path,filename]); % Data, time
[En] = wavefield_energy(Data);
E(1:length(En),c+1) = En;  
deltat= -3.9*1e-6;
%Delta=[0.0293,-0.02563,-0.02563];
for c=1:length(tasks)
    plot((t(:,c)+deltat)*1e3,E(:,c),lc(c),'LineWidth',1);
    hold on;
end
plot(time*1e3,En,'k','LineWidth',2);
axis([0 t_frames(end/2)*1e3 0 1]);

set(gcf,'color','white');
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
xlabel('t [ms]','Fontsize',12);
ylabel('E [-]','Fontsize',12);
hold on;
title(['Frequency ',num2str(f_2/1000),' [kHz]'],'Fontsize',12); 
legend({['$\alpha_{z}=$ ',num2str(alphaz(1)),', $\alpha_{xy}=$ ',num2str(alphaxy(1))],['$\alpha_{z}=$ ',num2str(alphaz(2)),', $\alpha_{xy}=$ ',num2str(alphaxy(2))],['$\alpha_{z}=$ ',num2str(alphaz(3)),', $\alpha_{xy}=$ ',num2str(alphaxy(3))],['$\alpha_{z}=$ ',num2str(alphaz(4)),', $\alpha_{xy}=$ ',num2str(alphaxy(4))],'experiment'},'Fontsize',12,'interpreter','latex');
figfilename = ['num_exp_energy',num2str(f_2/1000),'_kHz.png'];
print([figure_output_path,figfilename],'-dpng', '-r600'); 
close all;
%% 100 kHz
% numerical
tasks=[6,9,12,15];% 100 kHz
c=0;
alphaz=zeros(length(tasks),1);
alphaxy=zeros(length(tasks),1);
E=zeros(1024,length(tasks)+1)+NaN;
t=zeros(1024,length(tasks)+1);
for test_case = tasks
    c=c+1;
    input_no = input_file_no(test_case);
    run(fullfile('inputs',['input',num2str(input_no)]));
    alphaz(c) = etad_z;
    alphaxy(c) = etad_xy;
    interim_output_name = [model_interim_path,filesep,num2str(test_case),'_output',filesep];

      if(overwrite)
            load([interim_output_name,'flat_shell_Vz_',num2str(test_case),'_982x982bottom']);
            Data=6*Data(1:491,1:491,:);% Data
            load([interim_output_name,'t_frames']); % t_frames (time vector)

           [En] = wavefield_energy(Data);
           E(1:length(En),c) = En;      
           t(1:length(t_frames),c) = t_frames;
      else
          fprintf([modelname,' test case: %d already exist\n'], test_case);
      end
end

% experiment
filename = list{3}; %100 kHz
load([input_data_path,filename]); % Data, time
[En] = wavefield_energy(Data);
E(1:length(En),c+1) = En;  
deltat= 0;
%Delta=[0.0293,-0.02563,-0.02563];
for c=1:length(tasks)
    plot((t(:,c)+deltat)*1e3,E(:,c),lc(c),'LineWidth',1);
    hold on;
end
plot(time*1e3,En,'k','LineWidth',2);
axis([0 t_frames(end/2)*1e3 0 1]);

set(gcf,'color','white');
set(gca,'Fontsize',10,'linewidth',1);
set(gca,'FontName','Times');
xlabel('t [ms]','Fontsize',12);
ylabel('E [-]','Fontsize',12);
hold on;
title(['Frequency ',num2str(f_2/1000),' [kHz]'],'Fontsize',12); 
legend({['$\alpha_{z}=$ ',num2str(alphaz(1)),', $\alpha_{xy}=$ ',num2str(alphaxy(1))],['$\alpha_{z}=$ ',num2str(alphaz(2)),', $\alpha_{xy}=$ ',num2str(alphaxy(2))],['$\alpha_{z}=$ ',num2str(alphaz(3)),', $\alpha_{xy}=$ ',num2str(alphaxy(3))],['$\alpha_{z}=$ ',num2str(alphaz(4)),', $\alpha_{xy}=$ ',num2str(alphaxy(4))],'experiment'},'Fontsize',12,'interpreter','latex');
figfilename = ['num_exp_energy',num2str(f_2/1000),'_kHz.png'];
print([figure_output_path,figfilename],'-dpng', '-r600'); 