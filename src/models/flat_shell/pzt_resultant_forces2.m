function [Nfi,Mfi] = pzt_resultant_forces2(ce,theta,d31,d32,Voltage,h,hpzt)
% PZT_RESULTANT_FORCES2   Resultant membrane forces and bending moments at neutral plane of the plate
% 
% Syntax: [Nfi,Mfi] = pzt_resultant_forces2(ce,theta,d31,d32,Voltage,h,hpzt)
% 
% Inputs: 
%    ce - elastic stiffness matrix of pzt under constant electric field, dimensions [3,3], Units: N/m^2
%    theta - angle between principal direction of pzt material and global x axis, double, Units: deg
%    d31, d32 - stress/charge constants, double, Units: [C/N] or [m/V] 
% 
% Outputs: 
%    Nifi - Resultant membrane forces vector [Nx;Ny;Nxy], double, dimensions [3, 1], Units: N
%    Mfi - Resultant bending moments at neutral plane of the plate, [Mx;My;Mxy], double, dimensions [3, 1], Units: N*m
% 
% Example: 
%    [Nfi,Mfi] = pzt_resultant_forces2(ce,theta,d31,d32,Voltage,h,hpzt)
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

q11 = ce(1,1);
q12 = ce(1,2);
q22 = ce(2,2);
q66 = ce(6,6);

 c=cos(theta*pi/180);
 s=sin(theta*pi/180);
 
 
 Te = [c^2,    s^2,      c*s;
        s^2,    c^2,    -s*c;
        -2*c*s,   2*c*s,  c^2-s^2];

 dpzt=[d31;d32;0];
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=cos(theta*pi/180); n=sin(theta*pi/180);
% elastic properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s11=q11.*m.^4+2*(q12+2*q66).*m.^2.*n.^2+q22.*n.^4;
s12=(q11+q22-4*q66).*m.^2.*n.^2+q12.*(m.^4+n.^4);
s16=(q11-q12-2*q66).*m.^3.*n+(q12-q22+2*q66).*m.*n.^3;
s22=q11.*n.^4+2*(q12+2*q66).*m.^2.*n.^2+q22.*m.^4;
s26=(q11-q12-2*q66).*n.^3.*m+(q12-q22+2*q66).*n.*m.^3;
%s66=(q11+q22-2*q12-2*q66).*m.^2.*n.^2+q66.*(m.^4+n.^4);
s66=(q11+q22-2*q12).*m.^2.*n.^2+q66.*(m.^2-n.^2).^2;
Qhat=[s11,s12,s16;
          s12,s22,s26;
          s16,s26,s66];
Nfi=-Qhat*((Te*dpzt)*Voltage);
Mfi=Nfi*(h/2+hpzt/2);
%---------------------- END OF CODE---------------------- 

% ================ [pzt_resultant_forces2.m] ================  
