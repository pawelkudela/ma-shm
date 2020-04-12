function [Data] = wavefield_delta_t_correction(Data,time,deltat,A)
% WAVEFIELD_DELTA_T_CORRECTION   Removes the time difference deltat 
%    It uses 1D linear interpolation 
%    It is useful i.e. when experimental and numerical signals differs by constant deltat
%    Signal is also multiplied by A
% 
% Syntax: [Data] = wavefield_delta_t_correction(Data,time,deltat,A)
% 
% Inputs: 
%    Data - Wavefield data, double matrix, dimensions [nx,ny,nft]
%    nx - number of points in x direction
%    ny - number of points in y direction
%    nft - number of time steps
%    time - time vector, dimensions [nft, 1], Units: s 
%    deltat - time difference, double, Units: s 
%    A - constant multiplier, double
% 
% Outputs: 
%    Data - corrected data field 
% 
% Example: 
%    [Data] = wavefield_delta_t_correction(Data,time,deltat,A) 
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

Dataq=interp1(time,permute(Data,[3,2,1]),time-deltat); 
Data = A*permute(Dataq,[3,2,1]);

%---------------------- END OF CODE---------------------- 

% ================ [wavefield_delta_t_correction.m] ================  
