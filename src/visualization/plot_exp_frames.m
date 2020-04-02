function plot_exp_frames(Data,time,filename,selected_frames,figure_output_path,normalization,caxis_cut,ColorMapName,fig_width,fig_height)
% PLOT_EXP_FRAMES   plot selected frames of wavefield from experimental measurements
%    wavefield is defined on regular grid of points (meshgrid)
%    frames are saved as png images 
% 
% Syntax: plot_exp_frames(Data,selected_frames,figure_output_path,normalization,caxis_cut,ColorMapName,field_variable,motion,fig_width,fig_height)
% 
% Inputs: 
%    Data - Wavefield data on regular grid, double, dimensions [Nx, Ny, nFrames], Units: m, m/s or m/s^2
%             Nx - number of points in x direction
%             Ny - number of points in y direction
%             nFrames - number of frames
%    time - time vector [s]
%    filename - name of file containing Data
%    selected_frames - list of frame numbers to plot, integer 
%    figure_output_path - path to folder in which figures are stored
%    normalization - normalization to the highest value of the wavefield, logical
%    caxis_cut - coefficient of caxis for cutting the amplitude values displayed, 0 < double caxis_cut <=1
%    ColorMapName - colour map name for plotting, string (custom or standard matlab colormaps)
%    field_variable - string: 'displacement', 'velocity' or 'acceleration' 
%    fig_width - figure width in cm
%    fig_height - figure height in cm
% 
% 
% Example: 
% 
% Other m-files required: flat_shell_variable_names.m
% Subfunctions: none 
% MAT-files required: none 
% See also: 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

%    ColorMapName - map name: 'map_sunset', 'map_sunset_interp', 'map_burd', 'map_burd_interp'
%    'map_brewer', 'map_brewer_interp', 'map_iridescent', 'map_iridescent_interp', 'map_discrete_rainbow'
is_custom_map = strcmp(ColorMapName,'map_sunset') || strcmp(ColorMapName,'map_sunset_interp') || strcmp(ColorMapName,'map_burd') ...
    || strcmp(ColorMapName,'map_burd_interp') || strcmp(ColorMapName,'map_brewer') || strcmp(ColorMapName,'map_brewer_interp')...
    || strcmp(ColorMapName,'map_iridescent') || strcmp(ColorMapName,'map_iridescent_interp') || strcmp(ColorMapName,'map_discrete_rainbow') ; 
if(is_custom_map)
    [cmap] = custom_colour_maps(ColorMapName);
end

Smax=max(max(max(abs(Data))));
Smin=-Smax;
for n=selected_frames
    tstring = sprintf('%5.2f',time(n)*1e6); % time in micro secunds
    figfilename = [filename,'_','frame',num2str(n),'_',tstring,'.png'];
    if(~normalization)
        Smax=max(max(abs(Data(:,:,n))));
        Smin=-Smax;
    end
    figure(n);
    surf(squeeze(((Data(:,:,n)))));
    
    axis equal;view(2);
    set(gcf,'color','white');axis off;
    shading interp; 
    if(is_custom_map)
        colormap(cmap);
    else
        colormap(ColorMapName);
    end
    set(gcf,'Renderer','zbuffer');

    caxis([caxis_cut*Smin,caxis_cut*Smax]);
    set(gca, 'Position',[0 0 1. 1.]); % figure without axis and white border
    set(gcf, 'Units','centimeters', 'Position',[10 10 fig_width fig_height]); 
    % remove unnecessary white space
    %set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02));
    %fig.PaperPositionMode   = 'auto';
    set(gcf,'PaperPositionMode','auto');
        
    print([figure_output_path,figfilename],'-dpng', '-r600'); 
    close all;
end

%---------------------- END OF CODE---------------------- 

% ================ [plot_meshgrid_frames.m] ================  
