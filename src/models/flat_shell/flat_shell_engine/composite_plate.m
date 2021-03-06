function [rho, m11, J11, a11, a22, a12, a16, a26, a66, a44, a45, a55, b11,b12,b16,b22,b26,b66,d11, d12, d16, d22, d26, d66]=composite_plate(h, h1, h2, rhom, rhof, em, ef, nim, nif, vol, alpha, lay)
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
m11=0; J11=0;
a44=0; a45=0; a55=0; 
d11=0; d12=0; d16=0; d22=0; d26=0; d66=0; 
a11=0;a12=0;a22=0;a16=0;a26=0;a66=0;
b11=0;b12=0;b16=0;b22=0;b26=0;b66=0;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of composite.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%