function [I_G_dofs,I_L_dofs] = global_local_12baskets_5dofs_v2(I_G,I_L,mode)
% GLOBAL_LOCAL_12BASKETS_5DOFS_v2   convert nodal numbers into degrees of freedom numbers 
% 
% Syntax: [I_G_dofs,I_L_dofs] = global_local_12baskets_5dofs_v2(I_G,I_L,mode)
% 
% Inputs: 
%    I_G - global node numbers, integer, dimensions [m, 12]
%    I_L - local node numbers, interger, dimensions [m, 12]
%    mode - 'cpu' or 'gpu', string
% 
% Outputs: 
%    I_G_dofs - global degree of freedom numbers, integer, dimensions [5*m, 1]
%    I_L_dofs - local degree of freedom numbers, interger, dimensions [5*m, 1]
%         m = fen*NofElNodes/12
% Example: 
%    [I_G_dofs,I_L_dofs] = global_local_12baskets_5dofs_v2(I_G,I_L,mode)
% 
% Other m-files required: none 
% Subfunctions: none 
% MAT-files required: none 
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2 
% 

% Author: Pawel Kudela, D.Sc., Ph.D., Eng. 
% Institute of Fluid Flow Machinery Polish Academy of Sciences 
% Mechanics of Intelligent Structures Department 
% email address: pk@imp.gda.pl 
% Website: https://www.imp.gda.pl/en/research-centres/o4/o4z1/people/ 

%---------------------- BEGIN CODE---------------------- 

I_G_dofs = zeros(5*size(I_G,1),12);
I_G_dofs(1:5:end,:)=5*I_G-4;
I_G_dofs(2:5:end,:)=5*I_G-3;
I_G_dofs(3:5:end,:)=5*I_G-2;
I_G_dofs(4:5:end,:)=5*I_G-1;
I_G_dofs(5:5:end,:)=5*I_G;

I_L_dofs = zeros(5*size(I_L,1),12);
I_L_dofs(1:5:end,:)=5*I_L-4;
I_L_dofs(2:5:end,:)=5*I_L-3;
I_L_dofs(3:5:end,:)=5*I_L-2;
I_L_dofs(4:5:end,:)=5*I_L-1;
I_L_dofs(5:5:end,:)=5*I_L;

switch mode
    case 'gpu'
    I_G_dofs=gpuArray(I_G_dofs);   
    I_L_dofs=gpuArray(I_L_dofs);
    case 'cpu'
    
end

%---------------------- END OF CODE---------------------- 

% ================ [global_local_12baskets_5dofs_v2.m] ================  
