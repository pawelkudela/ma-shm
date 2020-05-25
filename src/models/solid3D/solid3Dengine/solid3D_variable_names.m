function [variable_name] = solid3D_variable_names(field_variable,motion)
% SOLID3D_VARIABLE_NAMES   Variable names specific to flat shell wavefield 
% 
% Syntax: [variable_name] = solid3D_variable_names(field_variable,motion) 
% 
% Inputs: 
%    field_variable - string: 'displacement', 'velocity' or 'acceleration' 
%    motion - integer defining motion type: 
%    1) Ux
%    2) Uy
%    3) Uz
%    8) sqrt((Ux).^2+(Uy).^2)
%    9) sqrt((Ux).^2+(Uy).^2 + Uz.^2)
% 
% Outputs: 
%    variable_name - name convention for saving, integer, dimensions [m, n], Units: - 
% 
% Example: 
%    [variable_name] = solid3D_variable_names(field_variable,motion) 
%    [variable_name] = solid3D_variable_names('velocity',3)  
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

switch field_variable
    case 'displacement'
        switch motion
            case 1
                variable_name = 'Ux';
            case 2
                variable_name = 'Uy';
            case 3
                variable_name = 'Uz';
            case 4
                variable_name = 'displacements_in_plane';
            case 5
                variable_name = 'total_displacements';
        end
    case 'velocity'
        switch motion
            case 1
                variable_name = 'Vx';
            case 2
                variable_name = 'Vy';
            case 3
                variable_name = 'Vz';
            case 4
                variable_name = 'velocities_in_plane';
            case 5
                variable_name = 'total_velocities';
        end
    case 'acceleration'
        switch motion
            case 1
                variable_name = 'Ax';
            case 2
                variable_name = 'Ay';
            case 3
                variable_name = 'Az';
            case 4
                variable_name = 'accelerations_in_plane';
            case 5
                variable_name = 'total_accelerations';
        end
end 

%---------------------- END OF CODE---------------------- 

% ================ [solid3D_variable_names.m] ================  
