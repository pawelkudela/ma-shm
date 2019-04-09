
Voltage = 10; % [V]
hpzt=0.1/1000; % pzt thickness [mm]
h=4/1000; % plate thickness [mm]
% pzt material properties Noliac NCE51
E1_pzt=11.5e10; %[N/m^2] youngs modulus at constant electric field
E2_pzt=E1_pzt; %[N/m^2] youngs modulus at constant electric field
ni12_pzt=0.32; % Poisson ratio  Noliac NCE51
ni21_pzt=0.32; % Poisson ratio  Noliac NCE51
G12_pzt =2.24e10; % shear modulus 
% Elastic stiffness matrix E
ce=[13.4, 8.89, 7.34,  0,   0,     0;
      8.89, 13.4, 7.34,  0,   0,     0;
      7.34, 7.34, 16.2,  0,   0,     0;
       0,     0,     0,    4.37, 0,    0;
       0,     0,     0,      0,   4.37, 0
       0,     0,     0,      0,   0,   2.24]*10e10; % [N/m^2]
% stress/charge constants
d31 = -208*10^-12; %[C/N], [m/V] Noliac NCE51
d32=d31;
%d32=-100*10^-12;
%theta - angle between principal direction of pzt material and global x axis, double, Units: deg
theta = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[Nfi,Mfi] = pzt_resultant_forces(ce,E1_pzt,E2_pzt,ni12_pzt,ni21_pzt,G12_pzt,theta,d31,d32,Voltage,h,hpzt)
[Nfi,Mfi] = pzt_resultant_forces2(ce,theta,d31,d32,Voltage,h,hpzt)