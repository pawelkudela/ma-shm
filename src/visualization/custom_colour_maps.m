function [cmap] = custom_colour_maps(ColorMapName)
% CUSTOM_COLOUR_MAPS   custom colour maps for surface plots 
% 
% Syntax: [cmap] = custom_colour_maps(ColorMapName)
% 
% Inputs: 
%    ColorMapName - map name: 'map_sunset', 'map_sunset_interp', 'map_burd', 'map_burd_interp'
%    'map_brewer', 'map_brewer_interp', 'map_iridescent', 'map_iridescent_interp', 'map_discrete_rainbow'
% 
% Outputs: 
%    cmap - colour map 3-column matrix
% 
% Example: 
%    [cmap] = custom_colour_maps(ColorMapName)
%    [cmap] = custom_colour_maps('sunset') 
% 
% Other m-files required: none 
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

N=256; % number of rows in colormap
x_interp = 0:1:N-1; %

switch ColorMapName
    % Diverging colour schemes
    % Diverging schemes are for ordered data between two extremes 
    % where the midpoint is important
    % colormaps that also works in colour-blind vision
    case 'map_sunset'
            cmap = [54,75,154;
                       74,123,183;
                       110,166,205;
                       152,202,225;
                       194,228,239;
                       234,236,204;
                       254,218,139;
                       253,179,102;
                       246,126,75;
                       221,61,45;
                       165,0,38];
    case 'map_burd'
          cmap = [33,102,172;
                     67,147,195;
                     146,197,222;
                     209,229,240;
                     247,247,247;
                     253,219,199;
                     244,165,130;
                     214,96,77;
                     178,24,43];
    case 'map_sunset_interp'
        map_sunset = [54,75,154;
                       74,123,183;
                       110,166,205;
                       152,202,225;
                       194,228,239;
                       234,236,204;
                       254,218,139;
                       253,179,102;
                       246,126,75;
                       221,61,45;
                       165,0,38];
        x1 = linspace(0,N-1,length(map_sunset));
        cmap = interp1(x1,map_sunset,x_interp','linear','extrap');
    case 'map_burd_interp'
        map_burd = [33,102,172;
                     67,147,195;
                     146,197,222;
                     209,229,240;
                     247,247,247;
                     253,219,199;
                     244,165,130;
                     214,96,77;
                     178,24,43];
        x2 = linspace(0,N-1,length(map_burd));
        cmap = interp1(x2,map_burd,x_interp','linear','extrap');
        
    % Sequential colour schemes
    % Sequential schemes are for ordered data from low to high
    case 'map_brewer'
          cmap = [255,255,229;
                     255,247,188;
                     254,227,145;
                     254,149,79;
                     251,154,41;
                     236,112,20;
                     204,76,2;
                     153,52,4;
                     102,37,6];
    case 'map_iridescent'
               cmap = [254,251,233;
                          254,247,213;
                          245,243,193;
                          234,240,181;
                          221,236,191;
                          208,231,202;
                          194,227,210;
                          181,221,216;
                          168,216,220;
                          155,210,225;
                          129,196,231;
                          123,188,231;
                          126,178,228;
                          136,165,221;
                          147,152,210;
                          155,138,196;
                          157,125,178;
                          154,112,158;
                          144,99,136;
                          128,87,112;
                          104,73,87;
                          70,53,58];
    case 'map_brewer_interp'
        map_brewer = [255,255,229;
                     255,247,188;
                     254,227,145;
                     254,149,79;
                     251,154,41;
                     236,112,20;
                     204,76,2;
                     153,52,4;
                     102,37,6];
        x3 = linspace(0,N-1,length(map_brewer));
        cmap = interp1(x3,map_brewer,x_interp','linear','extrap');
    case 'map_iridescent_interp'
        map_iridescent = [254,251,233;
                          254,247,213;
                          245,243,193;
                          234,240,181;
                          221,236,191;
                          208,231,202;
                          194,227,210;
                          181,221,216;
                          168,216,220;
                          155,210,225;
                          129,196,231;
                          123,188,231;
                          126,178,228;
                          136,165,221;
                          147,152,210;
                          155,138,196;
                          157,125,178;
                          154,112,158;
                          144,99,136;
                          128,87,112;
                          104,73,87;
                          70,53,58];
        x4 = linspace(0,N-1,length(map_iridescent));
        cmap = interp1(x4,map_iridescent,x_interp','linear','extrap');   
    case 'map_discrete_rainbow'
             cmap = [209,187,215;
                        174,118,163;
                        136,46,114;
                        25,101,176;
                        82,137,199;
                        123,175,222;
                        78,178,101;
                        144,201,135;
                        202,224,171;
                        247,240,86;
                        246,193,65;
                        241,147,45;
                        232,96,28;
                        220,5,12];
end
cmap = cmap/255;
%---------------------- END OF CODE---------------------- 

% ================ [custom_colour_maps.m] ================  
