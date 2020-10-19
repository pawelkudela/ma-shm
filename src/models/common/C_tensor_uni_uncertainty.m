% input
e11 = 125.5;
delta_e11 = 2.4;
ni12=0.37;
delata_ni12 = 0.08;
ni23=0.45;
delta_ni23=0.02;
e22=8.7;
delta_e22=0.1;
g12=4.135;
delta_g12 = 0.017*4.135; % 0.0703
%
e33=e22;
ni13=ni12;
g23=e22/(2*(1+ni23));
ni21=ni12*e22/e11;
ni31=ni13*e33/e11;
ni32=e33./e22.*ni23;
g13=g12;
Delta=(1-ni12.*ni21-ni13.*ni31-ni12.*ni23.*ni31-ni13.*ni21.*ni32-ni23.*ni32);
C11=e11.*(1-ni23.*ni32)./Delta;
C22=e22.*(1-ni31.*ni13)./Delta;
C33=e33.*(1-ni12.*ni21)./Delta;
C44=g23;
C55=g13;
C66=g12;
C12=(ni21+ ni31.*ni23).*e11./Delta;
C13=(ni31+ ni21.*ni32).*e11./Delta;
C23=(ni32+ ni12.*ni31).*e22./Delta;

% error propagation
delta_e11_percentage=delta_e11/e11;
delta_ni12_percentage=delata_ni12/ni12;
delta_ni13_percentage=delta_ni12_percentage;
delta_ni23_percentage=delta_ni23/ni23;
delta_e22_percentage=delta_e22/e22;
delta_e33_percentage=delta_e22_percentage;

delta_g23_percentage = sqrt(delta_e22_percentage^2+delta_ni23_percentage^2);
delta_g23 = g23*delta_g23_percentage;
delta_ni21_percentage = sqrt(delta_ni12_percentage^2+delta_e22_percentage^2+delta_e11_percentage^2);
delta_ni21 = ni21*delta_ni21_percentage;
delta_ni31 = delta_ni21;
delta_ni31_percentage = delta_ni31/ni31;
delta_ni32_percentage = sqrt(delta_e22_percentage^2 + delta_e22_percentage^2 + delta_ni23_percentage^2 );
delta_ni32 = ni32*delta_ni32_percentage;
delta_g13 = delta_g12;

%
% C11 130.00 \pm 12.22
% C12    6.06 \pm 1.08
% C13    6.06 \pm 1.08
% C22  11.19 \pm 3.51
% C23   5.19 \pm 0.42
% C33  11.19 \pm 3.51
% C44   3.00 \pm 0.14
% C55  4.13 \pm 0.07
% C66  4.13 \pm 0.07

delta_Delta = sqrt((ni12*ni21*sqrt(delata_ni12_percentage^2+delta_ni21_percentage^2))^2 +...
    (ni13*ni31*sqrt(delta_ni13_percentage^2 + delta_ni31_percentage^2))^2 +...
    (ni12.*ni23.*ni31*sqrt(delta_ni12_percentage^2 + delta_ni23_percentage^2) + delta_ni31_percentage^2)^2 +...
    (ni13.*ni21.*ni32*sqrt(delta_ni13_percentage^2 + delta_ni21_percentage^2) + delta_ni32_percentage^2)^2 +...
    (ni23.*ni32*sqrt(delta_ni23_percentage^2 + delta_ni32_percentage^2))^2);
delta_Delta_percentage = delta_Delta/Delta;

C11_delta_percentage = sqrt(delta_e11_percentage^2+delta_ni23_percentage^2+delta_ni32_percentage^2+delta_Delta_percentage^2);
C11_delta = C11*C11_delta_percentage;

C22_delta_percentage = sqrt(delta_e22_percentage^2+delta_ni31_percentage^2+delta_ni13_percentage^2+delta_Delta_percentage^2);
C22_delta = C22*C22_delta_percentage;

C33_delta_percentage = sqrt(delta_e33_percentage^2+delta_ni12_percentage^2+delta_ni21_percentage^2+delta_Delta_percentage^2);
C33_delta = C33*C33_delta_percentage;

C44_delta = delta_g23;
C55_delta = delta_g13;
C66_delta = delta_g12;

ni1= (ni21+ ni31.*ni23);
ni1_delta = sqrt(delta_ni21^2 + (ni31.*ni23*(sqrt(delta_ni31_percentage^2 + delta_ni23_percentage^2)))^2);
ni1_delta_percentage = ni1_delta/ni1; 
C12_delta_percentage = sqrt(ni1_delta_percentage^2+delta_e11_percentage^2+delta_Delta_percentage^2);
C12_delta = C12*C12_delta_percentage;

ni2=ni31+ ni21.*ni32;
ni2_delta = sqrt(delta_ni31^2 + (ni21.*ni32*(sqrt(delta_ni21_percentage^2 + delta_ni32_percentage^2)))^2);
ni2_delta_percentage = ni2_delta/ni2;
C13_delta_percentage = sqrt(ni2_delta_percentage^2+delta_e11_percentage^2+delta_Delta_percentage^2);
C13_delta = C13*C13_delta_percentage;

ni3=ni32+ ni12.*ni31;
ni3_delta = sqrt(delta_ni32^2 + (ni12.*ni31*(sqrt(delta_ni12_percentage^2 + delta_ni31_percentage^2)))^2);
ni3_delta_percentage = ni3_delta/ni3; 
C23_delta_percentage = sqrt(ni3_delta_percentage^2+delta_e22_percentage^2+delta_Delta_percentage^2);
C23_delta = C23*C23_delta_percentage;

