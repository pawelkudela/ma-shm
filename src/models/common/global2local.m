% matrix transforamtion global to local
function [l, m, n] = global2local(alpha, beta, gamma)
% alpha - x axis rotation; beta - y axis rotation; gamma - z axis rotation; 
    
    R_x = [1 0 0;
           0 cos(alpha) -sin(alpha);
           0 sin(alpha) cos(alpha)];
    R_y = [cos(beta) 0 sin(beta);
           0 1 0;
           -sin(beta) 0 cos(beta)];
    R_z = [cos(gamma) -sin(gamma) 0;
           sin(gamma) cos(gamma) 0;
           0 0 1];
                    
    R_xyz = round((R_x*R_y*R_z)*1e12)*1e-12;
    l = R_xyz(:,1);
    m = R_xyz(:,2);
    n = R_xyz(:,3);
end