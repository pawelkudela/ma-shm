function [Cmap] = colormapwlasny(whitethreshold,blackthreshold)
Cmap = 1-([blackthreshold:1/255:1-whitethreshold ; blackthreshold:1/255:1-whitethreshold ; blackthreshold:1/255:1-whitethreshold]');

