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
modelfolder = 'solid3D'; % name of folder

% prepare model input paths
model_riveted_delam90 = fullfile(projectroot,'src','models','solid3D','outputs','output1',filesep);
model_riveted_delam0 = fullfile(projectroot,'src','models','solid3D','outputs','output2',filesep);
model_no_rivets_delam90 = fullfile(projectroot,'src','models','solid3D','outputs','output3',filesep);
model_no_rivets_delam0 = fullfile(projectroot,'src','models','solid3D','outputs','output4',filesep);
model_baseline90 = fullfile(projectroot,'src','models','solid3D','outputs','output5',filesep);
model_baseline0 = fullfile(projectroot,'src','models','solid3D','outputs','output6',filesep);


% Input
nft=100000;    % total number of samples
tt=0.4e-3;   % total calculation time [s] % 
dt=tt/(nft);
t=dt:dt:tt;
t=t*1e3; % [ms]
% size 12cm by 8cm (1-column text)
fig_width = 12; fig_height = 10; 
linewidth = 0.5;

% load numerical data
    test_case = 1;
    load([model_riveted_delam90,'voltage',num2str(test_case)]);
    v_riveted_delam90=voltage;
    
    test_case = 2;
    load([model_riveted_delam0,'voltage',num2str(test_case)]);
    v_riveted_delam0=voltage;
    
    test_case = 3;
    load([model_no_rivets_delam90,'voltage',num2str(test_case)]);
    v_no_rivets_delam90=voltage;
    
    test_case = 4;
    load([model_no_rivets_delam0,'voltage',num2str(test_case)]);
    v_no_rivets_delam0=voltage;
    
    test_case = 5;
    load([model_baseline90,'voltage',num2str(test_case)]);
    v_baseline90=voltage;
    
    test_case = 6;
    load([model_baseline0,'voltage',num2str(test_case)]);
    v_baseline0=voltage;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% figures
 % effect of delamination
 % actuator no 3 'path_3_7';
 sensor=7;
  
    figname = 'path_3_7_delam';
    figfilename = fullfile(figure_output_path,figname);
    figure;subplot(2,1,1);
    set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height],'Color','w'); 
    hold on
    box on
    plot(t, v_riveted_delam90(:,sensor),'k:','LineWidth',4*linewidth);

    set(gca,'FontName','Times');
    hold on
    box on
    plot(t,v_baseline90(:,sensor),'r','LineWidth',2*linewidth);
    xlabel('t [ms]');
    ylabel('S [V]');
    legend('with delamination','no delamination','Location','NorthEast','Orientation','horizontal');
    %legend('differential exp.','differential num.','Location','NorthEastOutside');
    text(0.28,-0.6,{'Numerical - [90]_4'},'FontName','Times');
    title('Path 3-7');
    
    subplot(2,1,2);
    plot(t, v_riveted_delam0(:,sensor),'k:','LineWidth',4*linewidth);
    set(gca,'FontName','Times');
    hold on
    box on
    plot(t,v_baseline0(:,sensor),'r','LineWidth',2*linewidth);
    xlabel('t [ms]');
    ylabel('S [V]');
    legend('with delamination','no delamination','Location','NorthEast','Orientation','horizontal');
    %legend('differential exp.','differential num.','Location','NorthEastOutside');
    text(0.28,-0.3,{'Numerical - [0]_4'},'FontName','Times');
    title('Path 3-7');
    %
    set(gca,'FontName','Times');
    set(gcf,'PaperPositionMode','auto');
    print(figfilename,'-dpng', '-r600'); 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % effect of delamination
 % actuator no 3 'path_3_8';
 sensor=8;
  
    figname = 'path_3_8_delam';
    figfilename = fullfile(figure_output_path,figname);
    figure;subplot(2,1,1);
    set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height],'Color','w'); 
    hold on
    box on
    plot(t, v_riveted_delam90(:,sensor),'k:','LineWidth',4*linewidth);

    set(gca,'FontName','Times');
    hold on
    box on
    plot(t,v_baseline90(:,sensor),'r','LineWidth',2*linewidth);
    xlabel('t [ms]');
    ylabel('S [V]');
    legend('with delamination','no delamination','Location','NorthEast','Orientation','horizontal');
    %legend('differential exp.','differential num.','Location','NorthEastOutside');
    text(0.28,-0.5,{'Numerical - [90]_4'},'FontName','Times');
    title('Path 3-8');
    
    subplot(2,1,2);
    plot(t, v_riveted_delam0(:,sensor),'k:','LineWidth',4*linewidth);
    set(gca,'FontName','Times');
    hold on
    box on
    plot(t,v_baseline0(:,sensor),'r','LineWidth',2*linewidth);
    xlabel('t [ms]');
    ylabel('S [V]');
    legend('with delamination','no delamination','Location','NorthEast','Orientation','horizontal');
    %legend('differential exp.','differential num.','Location','NorthEastOutside');
    text(0.28,-0.1,{'Numerical - [0]_4'},'FontName','Times');
    title('Path 3-8');
    %
    set(gca,'FontName','Times');
    set(gcf,'PaperPositionMode','auto');
    print(figfilename,'-dpng', '-r600'); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% figures
 % effect of rivets
 % actuator no 3 'path_3_7';
 sensor=7;
  
    figname = 'path_3_7_rivets';
    figfilename = fullfile(figure_output_path,figname);
    figure;subplot(2,1,1);
    set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height],'Color','w'); 
    hold on
    box on
    plot(t, v_riveted_delam90(:,sensor),'k:','LineWidth',4*linewidth);

    set(gca,'FontName','Times');
    hold on
    box on
    plot(t,v_no_rivets_delam90(:,sensor),'r','LineWidth',2*linewidth);
    xlabel('t [ms]');
    ylabel('S [V]');
    legend('with rivets','no rivets','Location','NorthEast','Orientation','horizontal');
    %legend('differential exp.','differential num.','Location','NorthEastOutside');
    text(0.28,-1,{'Numerical - [90]_4'},'FontName','Times');
    title('Path 3-7');
    
    subplot(2,1,2);
    plot(t, v_riveted_delam0(:,sensor),'k:','LineWidth',4*linewidth);
    set(gca,'FontName','Times');
    hold on
    box on
    plot(t,v_no_rivets_delam0(:,sensor),'r','LineWidth',2*linewidth);
    xlabel('t [ms]');
    ylabel('S [V]');
    legend('with rivets','no rivets','Location','NorthEast','Orientation','horizontal');
    %legend('differential exp.','differential num.','Location','NorthEastOutside');
    text(0.28,-0.3,{'Numerical - [0]_4'},'FontName','Times');
    title('Path 3-7');
    %
    set(gca,'FontName','Times');
    set(gcf,'PaperPositionMode','auto');
    print(figfilename,'-dpng', '-r600'); 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % effect of rivets
 % actuator no 3 'path_3_8';
 sensor=8;
  
    figname = 'path_3_8_rivets';
    figfilename = fullfile(figure_output_path,figname);
    figure;subplot(2,1,1);
    set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height],'Color','w'); 
    hold on
    box on
    plot(t, v_riveted_delam90(:,sensor),'k:','LineWidth',4*linewidth);

    set(gca,'FontName','Times');
    hold on
    box on
    plot(t,v_no_rivets_delam90(:,sensor),'r','LineWidth',2*linewidth);
    xlabel('t [ms]');
    ylabel('S [V]');
    legend('with rivets','no rivets','Location','NorthEast','Orientation','horizontal');
    %legend('differential exp.','differential num.','Location','NorthEastOutside');
    text(0.28,-0.5,{'Numerical - [90]_4'},'FontName','Times');
    title('Path 3-8');
    
    subplot(2,1,2);
    plot(t, v_riveted_delam0(:,sensor),'k:','LineWidth',4*linewidth);
    set(gca,'FontName','Times');
    hold on
    box on
    plot(t,v_no_rivets_delam0(:,sensor),'r','LineWidth',2*linewidth);
    xlabel('t [ms]');
    ylabel('S [V]');
    legend('with rivets','no rivets','Location','NorthEast','Orientation','horizontal');
    %legend('differential exp.','differential num.','Location','NorthEastOutside');
    text(0.28,-0.1,{'Numerical - [0]_4'},'FontName','Times');
    title('Path 3-8');
    %
    set(gca,'FontName','Times');
    set(gcf,'PaperPositionMode','auto');
    print(figfilename,'-dpng', '-r600'); 