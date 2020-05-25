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
modelname_baseline =  'flat_shell_Jochen_signals_90_damping_baseline';
modelname_baseline2 =  'flat_shell_Jochen_signals_90_damping_baseline2';
modelname_added_mass = 'flat_shell_Jochen_signals_90_damping_added_mass';
modelname_added_mass2 = 'flat_shell_Jochen_signals_90_damping_added_mass2';
modelname_delamination = 'flat_shell_Jochen_signals_90_damping';
% prepare model input paths
model_baseline_path = prepare_model_paths('raw','num',modelfolder,modelname_baseline);
model_baseline_path2 = prepare_model_paths('raw','num',modelfolder,modelname_baseline2);
model_added_mass_path = prepare_model_paths('raw','num',modelfolder,modelname_added_mass);
model_added_mass_path2 = prepare_model_paths('raw','num',modelfolder,modelname_added_mass2);
model_delamination_path = prepare_model_paths('raw','num',modelfolder,modelname_delamination);
%% Input for flat_shell
tasks = [1:4]; % 
% scaling of numerical data
scaling=40;
%% input for experiment
% create path to the experimental data folder
input_data_path = fullfile( projectroot, 'data','external','exp','Jochen-Moll-Collab', 'SHM_plate',filesep );
% load filtered data
load([input_data_path,'data_all_filt']);
data_all=data_all_filt;
load([input_data_path,'data_all_baseline_filt']);
data_all_baseline = data_all_baseline_filt;
load([input_data_path,'time']);
time = time * 1e6;

% size 12cm by 8cm (1-column text)
fig_width = 12; fig_height = 10; 
linewidth = 0.5;
%% actuator no 1 'path_1_7';
actuator = 1; sensor = 7;
signal_no = 6; 
test_case = 1;
    % load numerical data
    model_path = [model_baseline_path,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    load([model_path,'time',num2str(test_case)]);
    v_baseline=voltage;
    t_baseline=t;
    model_path = [model_baseline_path2,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    v_baseline2=voltage;
    model_path = [model_added_mass_path,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    v_mass=voltage;
    model_path = [model_added_mass_path2,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    v_mass2=voltage;
    t_mass=t;
    model_path = [model_delamination_path,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    v_delam=voltage;
    t_delam=t;
    % 
 
    
    
figname = 'path_1_7';
figfilename = fullfile(figure_output_path,figname);
figure;

set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height],'Color','w'); % size 12cm by 8cm (1-column text)
subplot(2,1,1)
hold on
box on
plot(time, data_all_baseline(:,6),'LineWidth',2*linewidth);
plot(time,data_all(:,6),'r','LineWidth',linewidth);
text(20,0.08,'E/M crosstalk','FontName','Times');
rectangle('Position',[0,-0.05,100,0.1],'Curvature',[0.8,0.8],'LineStyle','--','LineWidth',1);
xlabel('Time [탎]');
ylabel('Amplitude [V]');
axis([0 1200 -0.17 0.17]);
%legend('baseline','added mass','Location','NorthEast');
legend('baseline', 'added mass','Location','NorthOutside','Orientation','horizontal');
title('Experimental transducer pairs 1-7');
set(gca,'FontName','Times');
subplot(2,1,2)
hold on
box on
plot(t_baseline*1e6,v_baseline(:,7)*scaling,'LineWidth',2*linewidth);
plot(t_delam*1e6,v_delam(:,7)*scaling,'g--','LineWidth',linewidth);
plot(t_mass*1e6,v_mass(:,7)*scaling,'r','LineWidth',linewidth);
xlabel('Time [탎]');
ylabel('Amplitude [V]');
axis([0 1200 -0.17 0.17])
%legend('baseline','delam', 'added mass','Location','NorthEast');
legend('baseline','delam', 'added mass','Location','North','Orientation','horizontal');
title('Numerical transducer pairs 1-7');
set(gca,'FontName','Times');
set(gcf,'PaperPositionMode','auto');

print(figfilename,'-dpng', '-r600'); 
%%
actuator = 2; sensor = 4;
signal_no = 13; 
test_case = 2;
    % load numerical data
    model_path = [model_baseline_path,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    load([model_path,'time',num2str(test_case)]);
    v_baseline=voltage;
    t_baseline=t;
    model_path = [model_baseline_path2,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    v_baseline2=voltage;
    model_path = [model_added_mass_path,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    v_mass=voltage;
    model_path = [model_added_mass_path2,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    v_mass2=voltage;
    t_mass=t;
    model_path = [model_delamination_path,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    v_delam=voltage;
    t_delam=t;
    
figname = 'path_2_4';
figfilename = fullfile(figure_output_path,figname);
figure;
set(gca,'FontName','Times');
set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height],'Color','w'); % size 12cm by 8cm (1-column text)
subplot(2,1,1)
hold on
box on
plot(time, data_all_baseline(:,13),'LineWidth',2*linewidth);
plot(time,data_all(:,13),'r','LineWidth',linewidth);
text(20,0.137,'E/M crosstalk','FontName','Times');
line([50,50],[0.05, 0.12],'LineWidth',1,'Color',[0 0 0]);
rectangle('Position',[0,-0.05,100,0.1],'Curvature',[0.8,0.8],'LineStyle','--','LineWidth',1);
xlabel('Time [탎]');
ylabel('Amplitude [V]');
axis([0 1200 -0.17 0.17]);
%legend('baseline','added mass','Location','NorthEast');
legend('baseline', 'added mass','Location','NorthOutside','Orientation','horizontal');
title('Experimental transducer pairs 2-4');
set(gca,'FontName','Times');
subplot(2,1,2)
hold on
box on
plot(t_baseline*1e6,v_baseline(:,4)*scaling,'LineWidth',2*linewidth);
plot(t_delam*1e6,v_delam(:,4)*scaling,'g--','LineWidth',linewidth);
plot(t_mass*1e6,v_mass(:,4)*scaling,'r','LineWidth',linewidth);
xlabel('Time [탎]');
ylabel('Amplitude [V]');
axis([0 1200 -0.17 0.17])
%legend('baseline','delam', 'added mass','Location','NorthEast');
legend('baseline','delam', 'added mass','Location','NorthEast','Orientation','horizontal');
title('Numerical transducer pairs 2-4');
set(gca,'FontName','Times');
set(gcf,'PaperPositionMode','auto');
print(figfilename,'-dpng', '-r600'); 

 
    %% actuator no 4 'path_4_8';
    test_case = 4;
    actuator = 4; sensor = 8;
    signal_no = 34; 
    % load numerical data
    model_path = [model_baseline_path,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    load([model_path,'time',num2str(test_case)]);
    v_baseline=voltage;
    t_baseline=t;
    model_path = [model_baseline_path2,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    v_baseline2=voltage;
    model_path = [model_added_mass_path,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    v_mass=voltage;
    model_path = [model_added_mass_path2,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    v_mass2=voltage;
    t_mass=t;
    model_path = [model_delamination_path,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    v_delam=voltage;
    t_delam=t;
    
    figname = 'path_4_8';
figfilename = fullfile(figure_output_path,figname);
figure;
set(gca,'FontName','Times');
set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height],'Color','w'); % size 12cm by 8cm (1-column text)
subplot(2,1,1)
hold on
box on
plot(time, data_all_baseline(:,34),'LineWidth',2*linewidth);
plot(time,data_all(:,34),'r','LineWidth',linewidth);
text(20,0.08,'E/M crosstalk','FontName','Times');
rectangle('Position',[0,-0.05,100,0.1],'Curvature',[0.8,0.8],'LineStyle','--','LineWidth',1);
xlabel('Time [탎]');
ylabel('Amplitude [V]');
axis([0 1200 -0.2 0.2]);
%legend('baseline','added mass','Location','NorthEast');
legend('baseline', 'added mass','Location','NorthOutside','Orientation','horizontal');
title('Experimental transducer pairs 4-8');
set(gca,'FontName','Times');
subplot(2,1,2)
hold on
box on
plot(t_baseline*1e6,v_baseline(:,8)*scaling,'LineWidth',2*linewidth);
plot(t_delam*1e6,v_delam(:,8)*scaling,'g--','LineWidth',linewidth);
plot(t_mass*1e6,v_mass(:,8)*scaling,'r','LineWidth',linewidth);
xlabel('Time [탎]');
ylabel('Amplitude [V]');
axis([0 1200 -0.2 0.2])
%legend('baseline','delam', 'added mass','Location','NorthEast');
legend('baseline','delam', 'added mass','Location','NorthEast','Orientation','horizontal');
title('Numerical transducer pairs 4-8');
set(gca,'FontName','Times');
set(gcf,'PaperPositionMode','auto');
print(figfilename,'-dpng', '-r600'); 
