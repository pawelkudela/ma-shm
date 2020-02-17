function [rho, m11, J11, a11, a22, a12, a16, a26, a66, a44, a45, a55, b11,b12,b16,b22,b26,b66,d11, d12, d16, d22, d26, d66]=composite_plate_pzt(h, h1, h2, rhom, rhof, em, ef, nim, nif, vol, alpha, lay,fen,NofElNodes,rho_pzt,pzt_thickness,Qpzt,pztEl)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mechanical properties of composite material components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gm=em./(1+nim)/2; gf=ef./(1+nif)/2;
rho=rhof.*vol+rhom.*(1-vol);
e11=ef.*vol+em.*(1-vol);
e22=((ef+em)+(ef-em).*vol)./((ef+em)-(ef-em).*vol).*em;

ni12=nif.*vol+nim.*(1-vol);

ni21=e22./e11.*ni12;

ni23=nif.*vol+nim.*(1-vol) .*(1.+nim-ni12.*em./e11) ./(1.-nim.^2 +nim.*ni12.*em./e11);

g12=((gf+gm)+(gf-gm).*vol)./((gf+gm)-(gf-gm).*vol).*gm;
g23=e22./2./(1+ni23);
g13=g12;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sierakowski page 46 eq. 2.33
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q11=e11./(1-ni12.*ni21);
q22=e22./(1-ni12.*ni21);
q44=g23;
q55=g13;
q66=g12;
q12=ni12.*e22./(1-ni12.*ni21);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global single layer properties 
% wykreslone sa: s13 s23 s33 s36, bo epsz=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphar=alpha*pi/180; m=cos(alphar); n=sin(alphar);
% elastic properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s11=q11.*m.^4+2*(q12+2*q66).*m.^2.*n.^2+q22.*n.^4;
s12=(q11+q22-4*q66).*m.^2.*n.^2+q12.*(m.^4+n.^4);
s16=(q11-q12-2*q66).*m.^3.*n+(q12-q22+2*q66).*m.*n.^3;
s22=q11.*n.^4+2*(q12+2*q66).*m.^2.*n.^2+q22.*m.^4;
s26=(q11-q12-2*q66).*n.^3.*m+(q12-q22+2*q66).*n.*m.^3;
s44=q44.*m.^2 + q55.*n.^2;
s45=(q55-q44).*m.*n;
s55=q55.*m.^2 + q44.*n.^2;
%s66=(q11+q22-2*q12-2*q66).*m.^2.*n.^2+q66.*(m.^4+n.^4);
s66=(q11+q22-2*q12).*m.^2.*n.^2+q66.*(m.^2-n.^2).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m11=zeros(fen*NofElNodes,1);
J11=zeros(fen*NofElNodes,1);
a44=zeros(fen*NofElNodes,1);
a45=zeros(fen*NofElNodes,1);
a55=zeros(fen*NofElNodes,1);
d11=zeros(fen*NofElNodes,1);
d12=zeros(fen*NofElNodes,1);
d16=zeros(fen*NofElNodes,1);
d22=zeros(fen*NofElNodes,1);
d26=zeros(fen*NofElNodes,1);
d66=zeros(fen*NofElNodes,1);
a11=zeros(fen*NofElNodes,1);
a12=zeros(fen*NofElNodes,1);
a22=zeros(fen*NofElNodes,1);
a16=zeros(fen*NofElNodes,1);
a26=zeros(fen*NofElNodes,1);
a66=zeros(fen*NofElNodes,1);
b11=zeros(fen*NofElNodes,1);
b12=zeros(fen*NofElNodes,1);
b16=zeros(fen*NofElNodes,1);
b22=zeros(fen*NofElNodes,1);
b26=zeros(fen*NofElNodes,1);
b66=zeros(fen*NofElNodes,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global multilayer properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:lay
% mass matrix components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  m11=m11+rho(i)*(h2(i)-h1(i));
  J11=J11+rho(i)*(h2(i)^3-h1(i)^3)/3;
% stiffness matrix components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  a44=a44+s44(i)*5/4*(h2(i)-h1(i)-4/3*(h2(i)^3-h1(i)^3)/h^2);
  a45=a45+s45(i)*5/4*(h2(i)-h1(i)-4/3*(h2(i)^3-h1(i)^3)/h^2);
  a55=a55+s55(i)*5/4*(h2(i)-h1(i)-4/3*(h2(i)^3-h1(i)^3)/h^2);
  
  a11=a11+s11(i)*(h2(i)-h1(i));
  a12=a12+s12(i)*(h2(i)-h1(i));
  a16=a16+s16(i)*(h2(i)-h1(i));
  a22=a22+s22(i)*(h2(i)-h1(i));
  a26=a26+s26(i)*(h2(i)-h1(i));
  a66=a66+s66(i)*(h2(i)-h1(i));
  
  b11=b11+s11(i)*(h2(i)^2-h1(i)^2)/2;
  b12=b12+s12(i)*(h2(i)^2-h1(i)^2)/2;
  b16=b16+s16(i)*(h2(i)^2-h1(i)^2)/2;
  b22=b22+s22(i)*(h2(i)^2-h1(i)^2)/2;
  b26=b26+s26(i)*(h2(i)^2-h1(i)^2)/2;
  b66=b66+s66(i)*(h2(i)^2-h1(i)^2)/2;
  
  d11=d11+s11(i)*(h2(i)^3-h1(i)^3)/3;
  d12=d12+s12(i)*(h2(i)^3-h1(i)^3)/3;
  d16=d16+s16(i)*(h2(i)^3-h1(i)^3)/3;
  d22=d22+s22(i)*(h2(i)^3-h1(i)^3)/3;
  d26=d26+s26(i)*(h2(i)^3-h1(i)^3)/3;
  d66=d66+s66(i)*(h2(i)^3-h1(i)^3)/3;
end
%
%% pzt layer 
pztEl=reshape(pztEl,1,[]);
ipztEl = zeros(length(pztEl)*NofElNodes,1);
c=0;
for ne=pztEl
    c=c+1;
    n1=(ne-1)*NofElNodes+1;
    n2=n1+NofElNodes-1;
    I=n1:n2;
    ipztEl((c-1)*NofElNodes+1:c*NofElNodes) = I; % indices of pzt nodes
end
% mass matrix components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m11(ipztEl)=m11(ipztEl)+rho_pzt*(pzt_thickness);
J11(ipztEl)=J11(ipztEl)+rho_pzt*((h/2+pzt_thickness/2)^3-(h/2)^3)/3;
 % stiffness matrix components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h2_=h/2+pzt_thickness/2;
h1_=h/2; 
a44(ipztEl)=a44(ipztEl)+Qpzt(4,4)*5/4*(h2_-h1_-4/3*(h2_^3-h1_^3)/(h+pzt_thickness)^2);
a45(ipztEl)=a45(ipztEl)+Qpzt(4,5)*5/4*(h2_-h1_-4/3*(h2_^3-h1_^3)/(h+pzt_thickness)^2);
a55(ipztEl)=a55(ipztEl)+Qpzt(5,5)*5/4*(h2_-h1_-4/3*(h2_^3-h1_^3)/(h+pzt_thickness)^2);

a11(ipztEl)=a11(ipztEl)+Qpzt(1,1)*(h2_-h1_);
a12(ipztEl)=a12(ipztEl)+Qpzt(1,2)*(h2_-h1_);
a16(ipztEl)=a16(ipztEl)+Qpzt(1,6)*(h2_-h1_);
a22(ipztEl)=a22(ipztEl)+Qpzt(2,2)*(h2_-h1_);
a26(ipztEl)=a26(ipztEl)+Qpzt(2,6)*(h2_-h1_);
a66(ipztEl)=a66(ipztEl)+Qpzt(6,6)*(h2_-h1_);

b11(ipztEl)=b11(ipztEl)+Qpzt(1,1)*(h2_^2-h1_^2)/2;
b12(ipztEl)=b12(ipztEl)+Qpzt(1,2)*(h2_^2-h1_^2)/2;
b16(ipztEl)=b16(ipztEl)+Qpzt(1,6)*(h2_^2-h1_^2)/2;
b22(ipztEl)=b22(ipztEl)+Qpzt(2,2)*(h2_^2-h1_^2)/2;
b26(ipztEl)=b26(ipztEl)+Qpzt(2,6)*(h2_^2-h1_^2)/2;
b66(ipztEl)=b66(ipztEl)+Qpzt(6,6)*(h2_^2-h1_^2)/2;

d11(ipztEl)=d11(ipztEl)+Qpzt(1,1)*(h2_^3-h1_^3)/3;
d12(ipztEl)=d12(ipztEl)+Qpzt(1,2)*(h2_^3-h1_^3)/3;
d16(ipztEl)=d16(ipztEl)+Qpzt(1,6)*(h2_^3-h1_^3)/3;
d22(ipztEl)=d22(ipztEl)+Qpzt(2,2)*(h2_^3-h1_^3)/3;
d26(ipztEl)=d26(ipztEl)+Qpzt(2,6)*(h2_^3-h1_^3)/3;
d66(ipztEl)=d66(ipztEl)+Qpzt(6,6)*(h2_^3-h1_^3)/3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of composite.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%