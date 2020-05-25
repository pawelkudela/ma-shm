% copy figures to paper folder and rename

load project_paths projectroot src_path;

figs_source_folder=[projectroot,'reports',filesep,'figures',filesep];
paper_folder = 'IOP-SMS';
fig_destination=[projectroot,'reports',filesep,'journal_papers',filesep,paper_folder,filesep];

modelname='Parallel_SMS_Num_Exp_Signals2'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figname=['path_1_7.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure8new.png'],'f');
figname=['path_2_4.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure9new.png'],'f');
figname=['path_4_8.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure10new.png'],'f');


modelname='Parallel_SMS_num_exp_energy_comparison'; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figname=['num_exp_energy16.5_kHz.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure12.png'],'f');
figname=['num_exp_energy50_kHz.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure13.png'],'f');
figname=['num_exp_energy100_kHz.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure14.png'],'f');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelname='Parallel_SMS_num_exp_signals_damp2'; 

figname=['num_exp_signals_damp_A_B_16.5_kHz.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure17.png'],'f');
figname=['num_exp_signals_damp_A_B_50_kHz.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure18.png'],'f');
figname=['num_exp_signals_damp_A_B_100_kHz.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure19.png'],'f');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelname='Parallel_SMS_Differential_Signals2';

figname=['path_1_7_diff.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure20new.png'],'f');
figname=['path_3_10_diff.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure21new.png'],'f');
figname=['path_4_8_diff.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure22new.png'],'f');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% experimental wavefield
modelname='Parallel_SMS_Jochen_wavefield_exp_colorbar';
figname=['colorbar_Vz_1_frame257_bottom.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure23a.png'],'f');

figname=['colorbar_Vz_1_frame401_bottom.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure23b.png'],'f');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% numerical wavefield
% delamination 
modelname='Parallel_SMS_Jochen_wavefield_90_colorbar';
figname=['colorbar_Vz_1_frame72_bottom.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,'1_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure23c.png'],'f');

figname=['colorbar_Vz_1_frame108_bottom.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,'1_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure23d.png'],'f');

% numerical wavefield
% added mass 
modelname='Parallel_SMS_Jochen_wavefield_90_added_mass_colorbar';
figname=['colorbar_Vz_1_frame72_bottom.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,'1_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure23e.png'],'f');

figname=['colorbar_Vz_1_frame108_bottom.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,'1_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure23f.png'],'f');

% numerical wavefield delamination vs added mass
modelname='Parallel_SMS_Jochen_wavefield_90_colorbar';
figname=['colorbar_cut_Vz_1_frame128_top.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,'1_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure24a.png'],'f');

modelname='Parallel_SMS_Jochen_wavefield_90_added_mass_colorbar';
figname=['colorbar_cut_Vz_1_frame128_top.png'];
fig_source=[figs_source_folder,modelname,'_out',filesep,'1_output',filesep,figname];
copyfile(fig_source,fig_destination);
movefile([fig_destination,figname],[fig_destination,'figure24b.png'],'f');
