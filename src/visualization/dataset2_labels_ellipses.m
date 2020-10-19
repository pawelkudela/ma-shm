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
modelfolder = pathstr(idx(end)+1:end); % name of folder
modelname = name; 
image_label_path1 = prepare_model_paths('raw','num','flat_shell','automesh_delam_rand'); % mesh parameters and labels
image_label_path2 = prepare_model_paths('raw','num','flat_shell','automesh_delam_rand2'); % mesh parameters and labels
image_label_path3 = prepare_model_paths('raw','num','flat_shell','automesh_delam_rand3'); % mesh parameters and labels
% prepare output paths
figure_output_path = prepare_figure_paths(modelname);
fig_width =5; % figure widht in cm
fig_height=5; % figure height in cm
figfilename = name;
%% prepare csv labels based on mesh parameters
% load mesh parameters
load([image_label_path1,filesep,'mesh_parameters']);
s=475;
position_no=zeros(s,1);
xCenter=zeros(s,1);
yCenter=zeros(s,1);
a=zeros(s,1);
b=zeros(s,1);
rotAngle=zeros(s,1);
for k=1:475%length(mesh_parameters)
    position_no(k) = k;
    xCenter(k) = mesh_parameters(k).xCenter;
    yCenter(k) = mesh_parameters(k).yCenter;
    a(k) = mesh_parameters(k).a;
    b(k) = mesh_parameters(k).b;
    rotAngle(k) = mesh_parameters(k).rotAngle;
end
load([image_label_path2,filesep,'mesh_parameters']);
for k=[257,318]
    position_no(k) = k;
    xCenter(k) = mesh_parameters(k).xCenter;
    yCenter(k) = mesh_parameters(k).yCenter;
    a(k) = mesh_parameters(k).a;
    b(k) = mesh_parameters(k).b;
    rotAngle(k) = mesh_parameters(k).rotAngle;
end
load([image_label_path3,filesep,'mesh_parameters']);
for k=[109,153,354]
    position_no(k) = k;
    xCenter(k) = mesh_parameters(k).xCenter;
    yCenter(k) = mesh_parameters(k).yCenter;
    a(k) = mesh_parameters(k).a;
    b(k) = mesh_parameters(k).b;
    rotAngle(k) = mesh_parameters(k).rotAngle;
end

%% plot ellipses
A=0.5;
for k=1:475
    
    xRadius=a(k);
    yRadius=b(k);
    plot_ellipse2(xCenter(k),yCenter(k),xRadius,yRadius,rotAngle(k),A);
    hold on;
end

set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
    set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
    % remove unnecessary white space
    %set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    %fig.PaperPositionMode   = 'auto';
    set(gcf,'PaperPositionMode','auto');
        
    print([figure_output_path,figfilename],'-dpng', '-r600'); 
    