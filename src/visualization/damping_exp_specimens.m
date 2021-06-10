% visualization of alpha damping parameter estimation for various specimens
% outputs from function: estimate_damping_coefficient

%    Algorithm is based on normalized energy of guided waves 
%    It uses simple grid search over damping coefficients alpha
%    and fits exponent function in the form f(t) = exp(-alpha*t) 
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
% size 12cm by 8cm (1-column text)
fig_width = 12; fig_height = 8; 
%% plain wave 50x50 cm specimens
L3_S4_B_freq=[16.5,50,50,75,75,100,100];% frequencies [kHz]
L3_S4_B_alpha=[1641,7831,7704,13223,13338,18798,18840];% alpha damping parameter

L3_S3_B_freq=[16.5,50,75,75,100,100];
L3_S3_B_alpha=[1579,7563,13091,12761,17691,17276];

L3_S2_B_freq=[16.5,50,100,100,100,100,150,200];
L3_S2_B_alpha=[1576,7428,17433,17426,17648,15067,21063,39272];

%% plain wave 120x120 quarter specimen
CFRP_120_plain_weave_quarter_freq=[16.5,50,100,150,200];
CFRP_120_plain_weave_quarter_alpha=[1724,11225,21804,42255,58590];
%% unidirectional 50x50 cm specimens
L1_S2_B_freq=[16.5,50,50,75,100,100];% frequencies [kHz]
L1_S2_B_alpha=[2518,9588,8591,13454,19132,18999];% alpha damping parameter

L1_S4_T_freq=[16.5,40,40,50,100];% frequencies [kHz]
L1_S4_T_alpha=[2373,7136,6856,9208,19541];% alpha damping parameter
%% unidirectional 120x120 quarter specimen
CFRP_120_UNI_quarter_freq=[16.5,50,100,150];
CFRP_120_UNI_quarter_alpha=[2368,9514,27210,45612];
%% plotting
% plain weave
figure;
plot(L3_S4_B_freq,L3_S4_B_alpha,'rx'); hold on;
plot(L3_S3_B_freq,L3_S3_B_alpha,'ro');
plot(L3_S2_B_freq,L3_S2_B_alpha,'rv');
plot(CFRP_120_plain_weave_quarter_freq,CFRP_120_plain_weave_quarter_alpha,'kd');

xlabel('f [kHz]','FontSize',12);
ylabel({'\alpha'},'FontSize',14);
xlim([0,220]);
legend('L3 S4 B','L3 S3 B','L3 S2 B','CFRP 120','Location','northwest');
title('CFRP plain weave');
set(gca,'FontName','Times');
set(gcf,'PaperPositionMode','auto');
figname = 'damping_plain_weave';
figfilename = fullfile(figure_output_path,figname);
print(figfilename,'-dpng', '-r600');

%
% unidirectional
figure;
plot(L1_S2_B_freq,L1_S2_B_alpha,'rx'); hold on;
plot(L1_S4_T_freq,L1_S4_T_alpha,'ro');

plot(CFRP_120_UNI_quarter_freq,CFRP_120_UNI_quarter_alpha,'kd');

xlabel('f [kHz]','FontSize',12);
ylabel({'\alpha'},'FontSize',14);
xlim([0,220]);
legend('L1 S2 B','L1 S4 T','CFRP 120 UNI','Location','northwest');
title('CFRP unidirectional');
set(gca,'FontName','Times');
set(gcf,'PaperPositionMode','auto');
figname = 'damping_unidirectional';
figfilename = fullfile(figure_output_path,figname);
print(figfilename,'-dpng', '-r600');

