function [Qpzt,epzt,gpzt]=pzt_const(Spzt,dp,epsT)

ep=dp*inv(Spzt);
gpzt=epsT-ep*dp';
%gpzt=epsT;
%elastic constants matrix of pzt element (stiffness matrix)
% change convention from Ex,Ey,Ez,Gxz,Gyz,Gxy into: Ex,Ey,Ez,Gxy,Gxz,Gyz
% E-epsilon, G-gamma 
Qpzt=[ Spzt(1,1) Spzt(1,2) Spzt(1,3) Spzt(1,6) Spzt(1,5)  Spzt(1,4)
       Spzt(2,1) Spzt(2,2) Spzt(2,3) Spzt(2,6) Spzt(2,5)  Spzt(2,4)
       Spzt(3,1) Spzt(3,2) Spzt(3,3) Spzt(3,6) Spzt(3,5)  Spzt(3,4)
       Spzt(6,1) Spzt(6,2) Spzt(6,3) Spzt(6,6) Spzt(6,5)  Spzt(6,4) 
       Spzt(5,1) Spzt(5,2) Spzt(5,3) Spzt(5,6) Spzt(5,5)  Spzt(5,4)
       Spzt(4,1) Spzt(4,2) Spzt(4,3) Spzt(4,6) Spzt(4,5)  Spzt(4,4)];
Spzt=Qpzt; % rearranged compliance matrix
Qpzt=inv(Spzt); % elastic constants matrix

% change convention from Ex,Ey,Ez,Gxz,Gyz,Gxy into: Ex,Ey,Ez,Gxy,Gxz,Gyz
epzt=[ ep(1,1) ep(1,2) ep(1,3) ep(1,6) ep(1,5) ep(1,4);
       ep(2,1) ep(2,2) ep(2,3) ep(2,6) ep(2,5) ep(2,4);
       ep(3,1) ep(3,2) ep(3,3) ep(3,6) ep(3,5) ep(3,4)];%[C/N]


