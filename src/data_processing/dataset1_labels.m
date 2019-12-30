clear all;close all;   warning off;clc;

load project_paths projectroot src_path;
%% Prepare output directories
% allow overwriting existing results if true
overwrite=false;
%overwrite=true;
% retrieve model name based on running file and folder
currentFile = mfilename('fullpath');
[pathstr,name,ext] = fileparts( currentFile );
idx = strfind( pathstr,filesep );
modelfolder = pathstr(idx(end)+1:end); % name of folder
modelname = name; 
image_label_path = prepare_model_paths('raw','num','flat_shell','automesh_delam1'); % mesh parameters and labels
% prepare output paths
dataset_output_path = prepare_data_processing_paths('processed','num',modelname);
%% prepare csv labels based on mesh parameters
% load mesh parameters
load([image_label_path,filesep,'mesh_parameters']);
position_no=zeros(length(mesh_parameters),1);
xCenter=zeros(length(mesh_parameters),1);
yCenter=zeros(length(mesh_parameters),1);
a=zeros(length(mesh_parameters),1);
b=zeros(length(mesh_parameters),1);
rotAngle=zeros(length(mesh_parameters),1);
for k=1:length(mesh_parameters)
    position_no(k) = mesh_parameters(k).position_no;
    xCenter(k) = mesh_parameters(k).xCenter;
    yCenter(k) = mesh_parameters(k).yCenter;
    a(k) = mesh_parameters(k).a;
    b(k) = mesh_parameters(k).b;
    rotAngle(k) = mesh_parameters(k).rotAngle;
end
position_no=cellstr(num2str(position_no,'%d'));
xCenter=cellstr(num2str(xCenter,'%.4f'));
yCenter=cellstr(num2str(yCenter,'%.4f'));
a=cellstr(num2str(a,'%.4f'));
b=cellstr(num2str(b,'%.4f'));
rotAngle=cellstr(num2str(rotAngle,'%.4f'));

T = table(position_no,xCenter,yCenter,a,b,rotAngle);
writetable(T,[dataset_output_path,filesep,'dataset1_labels.csv']);
%% move labels in the form of figures
status=movefile([image_label_path,filesep,'m_*'],[dataset_output_path,filesep]);
if (status)
    fprintf('Image label files has been moved succesfully\n');
else
    fprintf('Failed to move image labels\n');
end