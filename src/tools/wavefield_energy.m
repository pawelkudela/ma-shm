function [En] = wavefield_energy(Data)
%   WAVEFIELD_ENERGY   Calculates wavefield energy En(t) 
%    Energy is normalized to the max value 
% 
% Syntax: [En] = wavefield_energy(Data) 
% 
% Inputs: 
%    Data - Wavefield data, double matrix, dimensions [nx,ny,nft]
%    nx - number of points in x direction
%    ny - number of points in y direction
%    nft - number of time steps
% 
% Outputs: 
%    En - Energy, double, vector, dimensions [nft,1] Units: - 
% 
% Example: 
%    [En] = wavefield_energy(Data)
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

Data = Data.^2;
E= squeeze(sum(sum(Data,1),2));
En = E/max(E); % normalization

%---------------------- END OF CODE---------------------- 

% ================ [wavefield_energy.m] ================  
