% input
clc;
clear all; 
e11 = 125.5;
delta_e11 = 2.4;
ni12=0.37;
delta_ni12 = 0.08;
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
delta_ni12_percentage=delta_ni12/ni12;
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

delta_Delta = sqrt((ni12*ni21*sqrt(delta_ni12_percentage^2+delta_ni21_percentage^2))^2 +...
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

% homogenized properties
% CFRP plate
lay=16;
lh = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]*2/1000/lay;
alpha1 = [45,0,-45,90,-45,0,45,90,90,45,0,-45,90,-45,0,45];
rho=zeros(lay,1)+1571;
h=sum(lh);
h2=zeros(1,lay);
h1=zeros(1,lay);
h2(1)=h/2; h1(1)=h/2-lh(1);

for i=2:lay;
  h2(i)=h2(i-1)-lh(i-1);
  h1(i)=h1(i-1)-lh(i);
end;

[a11,a12,a13,a16,a22,a23,a33,a26,a36,a44,a45,a55,a66,J0]=composite3D_transf_C(h,h1,h2,rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,alpha1',lay);

% stringer

lay=12;
lh = [1,1,1,1,1,1,1,1,1,1,1,1]*1.5/1000/lay;
alpha2 = [-45,0,90,45,90,-45,-45,90,45,90,0,-45];
rho=zeros(lay,1)+1571;
h=sum(lh);
h2=zeros(1,lay);
h1=zeros(1,lay);
h2(1)=h/2; h1(1)=h/2-lh(1);

for i=2:lay;
  h2(i)=h2(i-1)-lh(i-1);
  h1(i)=h1(i-1)-lh(i);
end;


[d11,d12,d13,d16,d22,d23,d33,d26,d36,d44,d45,d55,d66,J0]=composite3D_transf_C(h,h1,h2,rho,C11,C12,C13,C22,C23,C33,C44,C55,C66,alpha2',lay);

% stiffener at angle 58 deg (rotation about x axis)
alpha=58*pi/180; beta = 0; gamma=0;
[l, m, n] = global2local(alpha, beta, gamma);
[D11,D12,D13,D14,D15,D16,D21,D22,D23,D24,D25,D26,D31,D32,D33,D34,D35,D36,D41,D42,D43,D44,D45,D46,D51,D52,D53,D54,D55,D56,D61,D62,D63,D64,D65,D66]=transform_D_Voigt(d11,d12,d13,d16,d22,d23,d26,d33,d36,d44,d45,d55,d66,l(1),l(2),l(3),m(1),m(2),m(3),n(1),n(2),n(3))
%
alpha=0;beta=0;
A11=zeros(lay,1);
for k=1:lay
    gamma=alpha2(k);
    [l, m, n] = global2local(alpha, beta, gamma);
    C16=0;C26=0;C36=0;C45=0;
    [D11,D12,D13,D14,D15,D16,D21,D22,D23,D24,D25,D26,D31,D32,D33,D34,D35,D36,D41,D42,D43,D44,D45,D46,D51,D52,D53,D54,D55,D56,D61,D62,D63,D64,D65,D66]=transform_D_Voigt(C11,C12,C13,C16,C22,C23,C26,C33,C36,C44,C45,C55,C66,l(1),l(2),l(3),m(1),m(2),m(3),n(1),n(2),n(3))
    A11(k)=D11;
end
for k=1:lay
    plot([A11(k),A11(k)],[h2(k),h1(k)]);
    hold on
end
[ksi,w]=gll(6);
plot(zeros(1,6)+d11,ksi*h/2,'ro');
%
B11=zeros(lay,1);
y=zeros(lay,1);
for k=1:lay
    B11(k)=A11(k);
    y(k)=(h2(k)+h1(k))/2;
end
 plot(B11,y,'g');
 zq=ksi(2:end-1)'*h/2;
 vq=interp1(y,B11,zq,'nearest');
 z=[-h/2;zq;h/2];
 v=[B11(1);vq;B11(lay)];
 hold on;
 plot(v,z,'md');