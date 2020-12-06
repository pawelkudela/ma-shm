function [MaxElementLength] = max_element_length(nodes,coords)
% MAX_ELEMENT_LENGTH   Find the longest edge in a 36-node mesh and output its length
% regular 6x6 node arrangement
% 
% Syntax: [MaxElementLength] = max_element_length(nodes,coords)
% 
% Inputs: 
%    nodes: incydency matrix of nodes; dimensions [fen,36]; fen - numeber of finite elements in a mesh
%    coords: coordinates, dimensions [m,2]; m - number of nodes in a mesh, Units: m 
% 
% Outputs: 
%    MaxElementLength - Length of the longest edge in the mesh, Units: mm
% 
% Example: 
%     [MaxElementLength] = max_element_length(nodes,coords)
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

MaxElementLength = max([max( abs(coords(nodes(:,1)) - coords(nodes(:,6)))),max( abs(coords(nodes(:,31)) - coords(nodes(:,36)))), max( abs(coords(nodes(:,6)) - coords(nodes(:,36)))),max( abs(coords(nodes(:,6)) - coords(nodes(:,36))))])*1e3;


%---------------------- END OF CODE---------------------- 

% ================ [max_element_length.m] ================  
