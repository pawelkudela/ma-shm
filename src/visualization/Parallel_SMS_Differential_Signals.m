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
modelname_baseline =  'flat_shell_Jochen_signals_90_damping_baseline2';
modelname_damaged = 'flat_shell_Jochen_signals_90_damping_added_mass2';
% prepare model input paths
model_baseline_path = prepare_model_paths('raw','num',modelfolder,modelname_baseline);
model_damaged_path = prepare_model_paths('raw','num',modelfolder,modelname_damaged);
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
fig_width = 12; fig_height = 4; 
linewidth = 0.5;
%% actuator no 1 'path_1_7_diff';
actuator = 1; sensor = 7;
signal_no = 6; 
test_case = 1;
    % load numerical data
    model_path = [model_baseline_path,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    load([model_path,'time',num2str(test_case)]);
    v_baseline=voltage;
    t_baseline=t;
    model_path = [model_damaged_path,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    v_mass=voltage;
    % differential signal
    figname = 'path_1_7_diff';
    figfilename = fullfile(figure_output_path,figname);
    figure;
    set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height],'Color','w'); % size 12cm by 4cm (1-column text)
    hold on
    box on
    plot(time, data_all_baseline(:,signal_no)-data_all(:,signal_no),'k:','LineWidth',4*linewidth);
    xlim([0 1200]);
    axis([0 1200 -0.01 0.01]);
    set(gca,'FontName','Times');
    hold on
    box on
    plot(t_baseline*1e6,v_baseline(:,7)*scaling-v_mass(:,7)*scaling,'r','LineWidth',2*linewidth);
    xlabel('t [µs]');
    ylabel('Amplitude [V]');
    xlim([0 1200]);
    axis([0 1200 -0.01 0.01])
    legend('differential exp.','differential num.','Location','NorthWest');
    title('Transducer pairs 1-7');
    set(gca,'FontName','Times');
    set(gcf,'PaperPositionMode','auto');
    print(figfilename,'-dpng', '-r600'); 
 %% actuator no 3 'path_3_10_diff';
    test_case = 3;
    actuator = 3; sensor = 10;
    signal_no = 28; 
    % load numerical data
    model_path = [model_baseline_path,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    load([model_path,'time',num2str(test_case)]);
    v_baseline=voltage;
    t_baseline=t;
    model_path = [model_damaged_path,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    v_mass=voltage;
    % differential signal
    figname = 'path_3_10_diff';
    figfilename = fullfile(figure_output_path,figname);
    figure;
    set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height],'Color','w'); % size 12cm by 4cm (1-column text)
    hold on
    box on
    plot(time, data_all_baseline(:,signal_no)-data_all(:,signal_no),'k:','LineWidth',4*linewidth);
    xlim([0 1200]);
    axis([0 1200 -0.01 0.01]);
    set(gca,'FontName','Times');
    hold on
    box on
    plot(t_baseline*1e6,v_baseline(:,10)*scaling-v_mass(:,10)*scaling,'r','LineWidth',2*linewidth);
    xlabel('t [µs]');
    ylabel('Amplitude [V]');
    xlim([0 1200]);
    axis([0 1200 -0.02 0.02])
    legend('differential exp.','differential num.','Location','NorthWest');
    title('Transducer pairs 3-10');
    set(gca,'FontName','Times');
    set(gcf,'PaperPositionMode','auto');
    print(figfilename,'-dpng', '-r600'); 
    %% actuator no 4 'path_4_8_diff';
    test_case = 4;
    actuator = 4; sensor = 8;
    signal_no = 34; 
    % load numerical data
    model_path = [model_baseline_path,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    load([model_path,'time',num2str(test_case)]);
    v_baseline=voltage;
    t_baseline=t;
    model_path = [model_damaged_path,filesep,num2str(test_case),'_output',filesep];
    load([model_path,'voltage',num2str(test_case)]);
    v_mass=voltage;
    % differential signal
    figname = 'path_4_8_diff';
    figfilename = fullfile(figure_output_path,figname);
    figure;
    set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height],'Color','w'); % size 12cm by 4cm (1-column text)
    hold on
    box on
    plot(time, data_all_baseline(:,signal_no)-data_all(:,signal_no),'k:','LineWidth',4*linewidth);
    xlim([0 1200]);
    axis([0 1200 -0.01 0.01]);
    set(gca,'FontName','Times');
    hold on
    box on
    plot(t_baseline*1e6,v_baseline(:,10)*scaling-v_mass(:,10)*scaling,'r','LineWidth',2*linewidth);
    xlabel('t [µs]');
    ylabel('Amplitude [V]');
    xlim([0 1200]);
    axis([0 1200 -0.02 0.02])
    legend('differential exp.','differential num.','Location','NorthWest');
    title('Transducer pairs 4-8');
    set(gca,'FontName','Times');
    set(gcf,'PaperPositionMode','auto');
    print(figfilename,'-dpng', '-r600'); 
