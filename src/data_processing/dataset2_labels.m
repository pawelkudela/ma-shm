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
dataset_output_path = prepare_data_processing_paths('processed','num',modelname);
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

position_no=cellstr(num2str(position_no,'%d'));
xCenter=cellstr(num2str(xCenter,'%.4f'));
yCenter=cellstr(num2str(yCenter,'%.4f'));
a=cellstr(num2str(a,'%.4f'));
b=cellstr(num2str(b,'%.4f'));
rotAngle=cellstr(num2str(rotAngle,'%.4f'));

T = table(position_no,xCenter,yCenter,a,b,rotAngle);
writetable(T,[dataset_output_path,filesep,'dataset2_labels.csv']);
%% move labels in the form of figures
status=movefile([image_label_path1,filesep,'*.png'],[dataset_output_path,filesep]);
if (status)
    fprintf('Image label files has been moved succesfully\n');
else
    fprintf('Failed to move image labels\n');
end

status=movefile([image_label_path2,filesep,'*.png'],[dataset_output_path,filesep]);
if (status)
    fprintf('Image label files has been moved succesfully\n');
else
    fprintf('Failed to move image labels\n');
end

status=movefile([image_label_path3,filesep,'*.png'],[dataset_output_path,filesep]);
if (status)
    fprintf('Image label files has been moved succesfully\n');
else
    fprintf('Failed to move image labels\n');
end