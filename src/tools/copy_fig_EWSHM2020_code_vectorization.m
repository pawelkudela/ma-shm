% copy figures to paper folder and rename

load project_paths projectroot src_path;

figs_source_folder=[projectroot,'reports',filesep,'figures',filesep];
paper_folder = 'EWSHM2020_code_vectorization';
fig_destination=[projectroot,'reports',filesep,'conference_papers',filesep,paper_folder,filesep,'figs',filesep];

modelname='flat_shell_EWSHM2020'; expname='EWSHM2020_exp_wavefield'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 16.5 kHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 1a num 16.5 kHz
figname=['Vz_1_frame64_bottom_249.74.png'];
fig_source=[figs_source_folder,'flat_shell',filesep,modelname,'_out',filesep,'1_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure1a.png'],'f');

% figure 1b exp 16.5 kHz
figname=['491x491p_16_5kHz_5HC_x5_15Vpp_frame65_250.00.png'];
fig_source=[figs_source_folder,expname,'_out',filesep,'1_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure1b.png'],'f');

% figure 1c num 16.5 kHz
figname=['Vz_1_frame128_bottom_499.50.png'];
fig_source=[figs_source_folder,'flat_shell',filesep,modelname,'_out',filesep,'1_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure1c.png'],'f');

% figure 1d exp 16.5 kHz
figname=['491x491p_16_5kHz_5HC_x5_15Vpp_frame129_500.00.png'];
fig_source=[figs_source_folder,expname,'_out',filesep,'1_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure1d.png'],'f');

% figure 1e num 16.5 kHz
figname=['Vz_1_frame192_bottom_749.25.png'];
fig_source=[figs_source_folder,'flat_shell',filesep,modelname,'_out',filesep,'1_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure1e.png'],'f');

% figure 1f exp 16.5 kHz
figname=['491x491p_16_5kHz_5HC_x5_15Vpp_frame193_750.00.png'];
fig_source=[figs_source_folder,expname,'_out',filesep,'1_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure1f.png'],'f');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 50 kHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 2a num 50 kHz
figname=['Vz_2_frame64_bottom_249.98.png'];
fig_source=[figs_source_folder,'flat_shell',filesep,modelname,'_out',filesep,'2_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure2a.png'],'f');

% figure 2b exp 50 kHz
figname=['492x492p_50kHz_5HC_x20_15Vpp_frame129_250.00.png'];
fig_source=[figs_source_folder,expname,'_out',filesep,'2_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure2b.png'],'f');

% figure 2c num 50 kHz
figname=['Vz_2_frame128_bottom_499.98.png'];
fig_source=[figs_source_folder,'flat_shell',filesep,modelname,'_out',filesep,'2_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure2c.png'],'f');

% figure 2d exp 50 kHz
figname=['492x492p_50kHz_5HC_x20_15Vpp_frame257_500.00.png'];
fig_source=[figs_source_folder,expname,'_out',filesep,'2_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure2d.png'],'f');

% figure 2e num 50 kHz
figname=['Vz_2_frame192_bottom_749.98.png'];
fig_source=[figs_source_folder,'flat_shell',filesep,modelname,'_out',filesep,'2_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure2e.png'],'f');

% figure 2f exp 50 kHz
figname=['492x492p_50kHz_5HC_x20_15Vpp_frame385_750.00.png'];
fig_source=[figs_source_folder,expname,'_out',filesep,'2_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure2f.png'],'f');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 100 kHz
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure 3a num 100 kHz
figname=['Vz_3_frame128_bottom_199.80.png'];
fig_source=[figs_source_folder,'flat_shell',filesep,modelname,'_out',filesep,'3_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure3a.png'],'f');

% figure 3b exp 100 kHz
figname=['491x491p_100kHz_5HC_x20_15Vpp_frame257_200.00.png'];
fig_source=[figs_source_folder,expname,'_out',filesep,'3_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure3b.png'],'f');

% figure 3c num 100 kHz
figname=['Vz_3_frame192_bottom_299.70.png'];
fig_source=[figs_source_folder,'flat_shell',filesep,modelname,'_out',filesep,'3_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure3c.png'],'f');

% figure 3d exp 100 kHz
figname=['491x491p_100kHz_5HC_x20_15Vpp_frame385_300.00.png'];
fig_source=[figs_source_folder,expname,'_out',filesep,'3_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure3d.png'],'f');

% figure 3e num 100 kHz
figname=['Vz_3_frame256_bottom_399.60.png'];
fig_source=[figs_source_folder,'flat_shell',filesep,modelname,'_out',filesep,'3_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure3e.png'],'f');

% figure 3f exp 100 kHz
figname=['491x491p_100kHz_5HC_x20_15Vpp_frame513_400.00.png'];
fig_source=[figs_source_folder,expname,'_out',filesep,'3_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure3f.png'],'f');