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
modelname =  'flat_shell_damping3';m=1;
% prepare model input/output paths
model_interim_path = prepare_model_paths('interim','num',modelfolder,modelname);

%% Input for flat_shell
input_file_no = [53:54]; % frequency = [50e3/5,100e3/5];
%% input for experiment
dataroot=fullfile('\\odroid-laser','laser','lamb_opt');
% create path to the experimental data folder
list = { [fullfile( dataroot, 'data','raw','exp', 'CFRP_120_UNI_quarter',filesep ),'492x492p_50kHz_5HC_x20_15Vpp'], ...
         [fullfile( dataroot, 'data','raw','exp', 'CFRP_120_UNI_quarter',filesep ),'492x492p_100kHz_5HC_x20_15Vpp'] ...
         [fullfile( dataroot, 'data','raw','exp', 'CFRP_120_quarter',filesep ),'497x497p_50kHz_5HC_x20_15Vpp'],...
         [fullfile( dataroot, 'data','raw','exp', 'CFRP_120_quarter',filesep ),'497x497p_100kHz_5HC_x20_15Vpp']};  
%load('Y:\lamb_opt\data\raw\exp\CFRP_120_UNI_quarter\492x492p_50kHz_5HC_x20_15Vpp');
%% input for figures
fig_width =5; % figure widht in cm
fig_height=5; % figure height in cm
lc='rgbm';


%% 50 kHz unidirectional
% numerical
tasks=[1];% 50 kHz
c=0;
alphaz=zeros(length(tasks),1);
alphaphi=zeros(length(tasks),1);
alphaxy=zeros(length(tasks),1);
E=zeros(1024,length(tasks)+1)+NaN;
t=zeros(1024,length(tasks)+1);
for test_case = tasks
    c=c+1;
    input_no = input_file_no(test_case);
    run(fullfile('inputs',['input',num2str(input_no)]));
    alphaz(c) = etad_z;
    alphaphi(c) = etad_z2;
    alphaxy(c) = etad_xy;
    interim_output_name = [model_interim_path,filesep,num2str(test_case),'_output',filesep];

      if(overwrite)
            load([interim_output_name,'flat_shell_Vz_',num2str(test_case),'_982x982bottom']);
            Data=Data(1:491,1:491,:);% Data
            load([interim_output_name,'t_frames']); % t_frames (time vector)

           [En] = wavefield_energy(Data);
           E(1:length(En),c) = En;   
           t(1:length(t_frames),c) = t_frames;
      else
          fprintf([modelname,' test case: %d already exist\n'], test_case);
      end
end

% experiment
filename = list{1}; %50 kHz
load(filename); % Data, time
[En] = wavefield_energy(Data);
E(1:length(En),c+1) = En;  
deltat= 0;

for c=1:length(tasks)
    if(c==3)   
        plot((t(:,c)+deltat)*1e3,E(:,c),'b--','LineWidth',1);
    else
        plot((t(:,c)+deltat)*1e3,E(:,c),lc(c),'LineWidth',1);
    end
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
if(m ==  2)
    legend({['$\alpha_{z}=$ ',num2str(alphaz(1))],...
        ['$\alpha_{z}=$ ',num2str(alphaz(2))],...
        ['$\alpha_{z}=$ ',num2str(alphaz(3))],...
        'experiment'},'Fontsize',12,'interpreter','latex');
    figfilename = ['num_exp_energy_uni2_',num2str(f_2/1000),'_kHz.png'];
else
 
%     legend({['$\alpha_{z}=$ ',num2str(alphaz(1)),', $\alpha_{xy}=$ ',num2str(alphaxy(1)),', $\alpha_{\phi}=$ ',num2str(alphaphi(1))],...
%         ['$\alpha_{z}=$ ',num2str(alphaz(2)),', $\alpha_{xy}=$ ',num2str(alphaxy(2)),', $\alpha_{\phi}=$ ',num2str(alphaphi(2))],...
%         ['$\alpha_{z}=$ ',num2str(alphaz(3)),', $\alpha_{xy}=$ ',num2str(alphaxy(3)),', $\alpha_{\phi}=$ ',num2str(alphaphi(3))],...
%         'experiment'},'Fontsize',12,'interpreter','latex');
legend({['$\alpha_{z}=$ ',num2str(alphaz(1)),', $\alpha_{xy}=$ ',num2str(alphaxy(1)),', $\alpha_{\phi}=$ ',num2str(alphaphi(1))],...
        'experiment'},'Fontsize',12,'interpreter','latex');
    figfilename = ['num_exp_energy_uni',num2str(f_2/1000),'_kHz.png'];
end
print([figure_output_path,figfilename],'-dpng', '-r600'); 
close all;
% errors calculation
% L2 error norm
L2er=zeros(length(tasks),1);
for c=1:length(tasks)
    L2er(c)=sqrt(sum((En(1:2:2*length(t_frames))-E(1:length(t_frames),c)).^2));   
end

% L2 relative error norm
L2er_relative=zeros(length(tasks),1);
for c=1:length(tasks)
    L2er_relative(c)=sqrt(sum((En(1:2:2*length(t_frames))-E(1:length(t_frames),c)).^2)/sum(En(1:length(t_frames)).^2));   
end
[L2er_relative]
% v1
% 0.1686
% 0.2061
% 0.1430
%v2
%  0.1545
%  0.1879
%  0.1341
%% 100 kHz unidirectional
% numerical
tasks=[2];% 100 kHz
c=0;
alphaz=zeros(length(tasks),1);
alphaphi=zeros(length(tasks),1);
alphaxy=zeros(length(tasks),1);
E=zeros(1024,length(tasks)+1)+NaN;
t=zeros(1024,length(tasks)+1);
for test_case = tasks
    c=c+1;
    input_no = input_file_no(test_case);
    run(fullfile('inputs',['input',num2str(input_no)]));
    alphaz(c) = etad_z;
    alphaxy(c) = etad_xy;
    alphaphi(c) = etad_z2;
    interim_output_name = [model_interim_path,filesep,num2str(test_case),'_output',filesep];

      if(overwrite)
            load([interim_output_name,'flat_shell_Vz_',num2str(test_case),'_982x982bottom']);
            Data=Data(1:491,1:491,:);% Data
            load([interim_output_name,'t_frames']); % t_frames (time vector)

           [En] = wavefield_energy(Data);
           E(1:length(En),c) = En;      
           t(1:length(t_frames),c) = t_frames;
      else
          fprintf([modelname,' test case: %d already exist\n'], test_case);
      end
end

% experiment
filename = list{2}; %100 kHz
load(filename); % Data, time
[En] = wavefield_energy(Data);
E(1:length(En),c+1) = En;  
deltat= 0;
%Delta=[0.0293,-0.02563,-0.02563];
%a=0.8;
a=1;
for c=1:length(tasks)
    if(c==3)
        plot((t(:,c)+deltat)*1e3,a*E(:,c),'b--','LineWidth',1);
    else
        plot((t(:,c)+deltat)*1e3,a*E(:,c),lc(c),'LineWidth',1);     
    end
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
if(m ==  2)
    legend({['$\alpha_{z}=$ ',num2str(alphaz(1))],...
        ['$\alpha_{z}=$ ',num2str(alphaz(2))],...
        ['$\alpha_{z}=$ ',num2str(alphaz(3))],...
        'experiment'},'Fontsize',12,'interpreter','latex');
    figfilename = ['num_exp_energy_uni2_',num2str(f_2/1000),'_kHz.png'];
else
%     legend({['$\alpha_{z}=$ ',num2str(alphaz(1)),', $\alpha_{xy}=$ ',num2str(alphaxy(1)),', $\alpha_{\phi}=$ ',num2str(alphaphi(1))],...
%         ['$\alpha_{z}=$ ',num2str(alphaz(2)),', $\alpha_{xy}=$ ',num2str(alphaxy(2)),', $\alpha_{\phi}=$ ',num2str(alphaphi(2))],...
%         ['$\alpha_{z}=$ ',num2str(alphaz(3)),', $\alpha_{xy}=$ ',num2str(alphaxy(3)),', $\alpha_{\phi}=$ ',num2str(alphaphi(3))],...
%         'experiment'},'Fontsize',12,'interpreter','latex');
    legend({['$\alpha_{z}=$ ',num2str(alphaz(1)),', $\alpha_{xy}=$ ',num2str(alphaxy(1)),', $\alpha_{\phi}=$ ',num2str(alphaphi(1))],...
        'experiment'},'Fontsize',12,'interpreter','latex');
    figfilename = ['num_exp_energy_uni',num2str(f_2/1000),'_kHz.png'];
end
print([figure_output_path,figfilename],'-dpng', '-r600'); 

% errors calculation
% L2 error norm
L2er=zeros(length(tasks),1);
for c=1:length(tasks)
    L2er(c)=sqrt(sum((En(1:2:2*length(t_frames))-E(1:length(t_frames),c)).^2));   
end

% L2 relative error norm
L2er_relative=zeros(length(tasks),1);
for c=1:length(tasks)
    L2er_relative(c)=sqrt(sum((En(1:2:2*length(t_frames))-E(1:length(t_frames),c)).^2)/sum(En(1:length(t_frames)).^2));   
end
[L2er_relative]
%v1
%  0.1686
%  0.2061
%  0.1430
%v2
%  0.1184
%  0.1378
%  0.1170


