function Cnew=damage_model(C,damage_degree)
% Apply damage degree to C tensor
% C - tensor of elastic constants (6x6)
% damage_degree - double in range 0:1 (0 no damage)
% Cnew - tensor of elastic constants with incorporated damage degree



% damage_degree = 0.5;
% C=[52.55 6.51    5.94     0       0     0;
%     6.51   51.83   5.88     0       0     0; 
%     5.94   5.88     10.28   0       0     0; 
%     0        0          0       2.93   0     0;
%     0        0          0       0       2.92 0;
%     0        0          0       0       0   3.81];
% Cinv=inv(C);
% Cd11 = damage_degree.*((1/Cinv(3,3))*(C(5,5)/C(3,3))^2 - (C(5,5)/C(3,3))^2);
% Cd22 = damage_degree.*((1/Cinv(3,3))*(C(4,4)/C(3,3))^2 - (C(4,4)/C(3,3))^2);
Cd11 = damage_degree.*((C(3,3))*(C(5,5)/C(3,3))^2 - (C(5,5)/C(3,3))^2);
Cd22 = damage_degree.*((C(3,3))*(C(4,4)/C(3,3))^2 - (C(4,4)/C(3,3))^2);
Cd33 = damage_degree.*C(3,3);
Cd44 = damage_degree.*C(4,4);
Cd55 = damage_degree.*C(5,5);
%Cd66 = damage_degree.*( ( 1/Cinv(3,3) )*(C(4,4)/C(3,3))*(C(5,5)/C(3,3)) - C(4,4)*C(5,5)/C(3,3)^2 );
Cd66 = damage_degree.*( C(3,3)*( C(4,4)/C(3,3) )*( C(5,5)/C(3,3) ) - (C(4,4)*C(5,5)/C(3,3)^2) );
Cd32 = damage_degree.*C(3,2);
Cd31 = damage_degree.*C(3,1);
Cd21 = (Cd22-Cd66)/2;

% figure
% plot(damage_degree,C(1,1)-Cd11);
% figure
% plot(damage_degree,C(2,2)-Cd22);
% figure
% plot(damage_degree,C(6,6)-Cd66);
% figure
% plot(damage_degree,C(2,1)-Cd21);


Cd = [Cd11 Cd21 Cd31   0       0      0;
        Cd21 Cd22 Cd32   0       0      0;
        Cd31 Cd32 Cd33   0       0      0;
        0         0      0   Cd44    0      0;
        0         0      0     0    Cd55    0;
        0         0      0     0      0    Cd66];
Cnew=C-Cd;

