clear all;close all;   warning off;clc;

load project_paths projectroot src_path;

% input
input_data_path = fullfile( projectroot, 'data','processed','num', filesep );
folder = 'Adaptive_filtering_wavefield_dataset2_out';
input = [input_data_path,folder,filesep,'IoU.mat'];
load(input);

% output 
output_data_path = fullfile( projectroot, 'data','processed','num', filesep );
folder = 'Adaptive_filtering_wavefield_dataset2_out';
output = [output_data_path,folder,filesep,'Adaptive_filtering_IoU.csv'];

IoU_mean_top = mean(IoU(:,1))
IoU_mean_bottom = mean(IoU(:,2))

IoU_best_top = max(IoU(:,1))
IoU_best_bottom = max(IoU(:,2))

IoU_worst_top = min(IoU(:,1))
IoU_worst_bottom = min(IoU(:,2))

No_of_undetected_top = 475-length(find(IoU(:,1)))
No_of_undetected_bottom = 475-length(find(IoU(:,2)))

% undetected case numbers at the bottom
I=find(IoU(:,2)==0)

T = table(IoU_worst_top, IoU_worst_bottom, IoU_best_top, IoU_best_bottom, IoU_mean_top,IoU_mean_bottom, No_of_undetected_top, No_of_undetected_bottom );
writetable(T,output);