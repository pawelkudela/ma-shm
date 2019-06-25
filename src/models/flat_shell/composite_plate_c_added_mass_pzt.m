function [rho, m11, J11, a11, a22, a12, a16, a26, a66, a44, a45, a55, b11,b12,b16,b22,b26,b66,d11, d12, d16, d22, d26, d66]=composite_plate_c_added_mass_pzt(h, h1, h2, rho, q11,q12,q22,q44,q55,q66,alpha,lay,delamEl,den_above,den_under,delamination_layer,fen,NofElNodes,rho_pzt,pzt_thickness,Qpzt,pztEl,addedMassEl,rho_added_mass,added_mass_thickness,Qm)


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
%% undelaminated region
unDelamEl=setdiff(1:fen,delamEl);
iunDelamEl = zeros(length(unDelamEl)*NofElNodes,1);
c=0;
for ne=unDelamEl
    c=c+1;
    n1=(ne-1)*NofElNodes+1;
    n2=n1+NofElNodes-1;
    I=n1:n2;
    iunDelamEl((c-1)*NofElNodes+1:c*NofElNodes) = I; % indices of undelaminated nodes
end
for i=1:lay
% mass matrix components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  m11(iunDelamEl)=m11(iunDelamEl)+rho(i)*(h2(i)-h1(i));
  J11(iunDelamEl)=J11(iunDelamEl)+rho(i)*(h2(i)^3-h1(i)^3)/3;
% stiffness matrix components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  a44(iunDelamEl)=a44(iunDelamEl)+s44(i)*5/4*(h2(i)-h1(i)-4/3*(h2(i)^3-h1(i)^3)/h^2);
  a45(iunDelamEl)=a45(iunDelamEl)+s45(i)*5/4*(h2(i)-h1(i)-4/3*(h2(i)^3-h1(i)^3)/h^2);
  a55(iunDelamEl)=a55(iunDelamEl)+s55(i)*5/4*(h2(i)-h1(i)-4/3*(h2(i)^3-h1(i)^3)/h^2);
  
  a11(iunDelamEl)=a11(iunDelamEl)+s11(i)*(h2(i)-h1(i));
  a12(iunDelamEl)=a12(iunDelamEl)+s12(i)*(h2(i)-h1(i));
  a16(iunDelamEl)=a16(iunDelamEl)+s16(i)*(h2(i)-h1(i));
  a22(iunDelamEl)=a22(iunDelamEl)+s22(i)*(h2(i)-h1(i));
  a26(iunDelamEl)=a26(iunDelamEl)+s26(i)*(h2(i)-h1(i));
  a66(iunDelamEl)=a66(iunDelamEl)+s66(i)*(h2(i)-h1(i));
  
  b11(iunDelamEl)=b11(iunDelamEl)+s11(i)*(h2(i)^2-h1(i)^2)/2;
  b12(iunDelamEl)=b12(iunDelamEl)+s12(i)*(h2(i)^2-h1(i)^2)/2;
  b16(iunDelamEl)=b16(iunDelamEl)+s16(i)*(h2(i)^2-h1(i)^2)/2;
  b22(iunDelamEl)=b22(iunDelamEl)+s22(i)*(h2(i)^2-h1(i)^2)/2;
  b26(iunDelamEl)=b26(iunDelamEl)+s26(i)*(h2(i)^2-h1(i)^2)/2;
  b66(iunDelamEl)=b66(iunDelamEl)+s66(i)*(h2(i)^2-h1(i)^2)/2;
  
  d11(iunDelamEl)=d11(iunDelamEl)+s11(i)*(h2(i)^3-h1(i)^3)/3;
  d12(iunDelamEl)=d12(iunDelamEl)+s12(i)*(h2(i)^3-h1(i)^3)/3;
  d16(iunDelamEl)=d16(iunDelamEl)+s16(i)*(h2(i)^3-h1(i)^3)/3;
  d22(iunDelamEl)=d22(iunDelamEl)+s22(i)*(h2(i)^3-h1(i)^3)/3;
  d26(iunDelamEl)=d26(iunDelamEl)+s26(i)*(h2(i)^3-h1(i)^3)/3;
  d66(iunDelamEl)=d66(iunDelamEl)+s66(i)*(h2(i)^3-h1(i)^3)/3;
end
%% elements above delamination
for j=1:length(den_above)
    DelamEl=den_above{j}; % list of elements above delamination no j
    iDelamEl = zeros(length(DelamEl)*NofElNodes,1);
    c=0;
    for ne=DelamEl
        c=c+1;
        n1=(ne-1)*NofElNodes+1;
        n2=n1+NofElNodes-1;
        I=n1:n2;
        iDelamEl((c-1)*NofElNodes+1:c*NofElNodes) = I; % indices of undelaminated nodes
    end
    for i=1:delamination_layer(j) % layers above delamination
    % mass matrix components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      m11(iDelamEl)=m11(iDelamEl)+rho(i)*(h2(i)-h1(i));
      J11(iDelamEl)=J11(iDelamEl)+rho(i)*(h2(i)^3-h1(i)^3)/3;
    % stiffness matrix components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      a44(iDelamEl)=a44(iDelamEl)+s44(i)*5/4*(h2(i)-h1(i)-4/3*(h2(i)^3-h1(i)^3)/h^2);
      a45(iDelamEl)=a45(iDelamEl)+s45(i)*5/4*(h2(i)-h1(i)-4/3*(h2(i)^3-h1(i)^3)/h^2);
      a55(iDelamEl)=a55(iDelamEl)+s55(i)*5/4*(h2(i)-h1(i)-4/3*(h2(i)^3-h1(i)^3)/h^2);

      a11(iDelamEl)=a11(iDelamEl)+s11(i)*(h2(i)-h1(i));
      a12(iDelamEl)=a12(iDelamEl)+s12(i)*(h2(i)-h1(i));
      a16(iDelamEl)=a16(iDelamEl)+s16(i)*(h2(i)-h1(i));
      a22(iDelamEl)=a22(iDelamEl)+s22(i)*(h2(i)-h1(i));
      a26(iDelamEl)=a26(iDelamEl)+s26(i)*(h2(i)-h1(i));
      a66(iDelamEl)=a66(iDelamEl)+s66(i)*(h2(i)-h1(i));

      b11(iDelamEl)=b11(iDelamEl)+s11(i)*(h2(i)^2-h1(i)^2)/2;
      b12(iDelamEl)=b12(iDelamEl)+s12(i)*(h2(i)^2-h1(i)^2)/2;
      b16(iDelamEl)=b16(iDelamEl)+s16(i)*(h2(i)^2-h1(i)^2)/2;
      b22(iDelamEl)=b22(iDelamEl)+s22(i)*(h2(i)^2-h1(i)^2)/2;
      b26(iDelamEl)=b26(iDelamEl)+s26(i)*(h2(i)^2-h1(i)^2)/2;
      b66(iDelamEl)=b66(iDelamEl)+s66(i)*(h2(i)^2-h1(i)^2)/2;

      d11(iDelamEl)=d11(iDelamEl)+s11(i)*(h2(i)^3-h1(i)^3)/3;
      d12(iDelamEl)=d12(iDelamEl)+s12(i)*(h2(i)^3-h1(i)^3)/3;
      d16(iDelamEl)=d16(iDelamEl)+s16(i)*(h2(i)^3-h1(i)^3)/3;
      d22(iDelamEl)=d22(iDelamEl)+s22(i)*(h2(i)^3-h1(i)^3)/3;
      d26(iDelamEl)=d26(iDelamEl)+s26(i)*(h2(i)^3-h1(i)^3)/3;
      d66(iDelamEl)=d66(iDelamEl)+s66(i)*(h2(i)^3-h1(i)^3)/3;
    end
end
%% elements under delamination
for j=1:length(den_under)
    DelamEl=den_under{j}; % list of elements under delamination no j
    iDelamEl = zeros(length(DelamEl)*NofElNodes,1);
    c=0;
    for ne=DelamEl
        c=c+1;
        n1=(ne-1)*NofElNodes+1;
        n2=n1+NofElNodes-1;
        I=n1:n2;
        iDelamEl((c-1)*NofElNodes+1:c*NofElNodes) = I; % indices of undelaminated nodes
    end
    for i= delamination_layer(j)+1:lay % layers under delamination
    % mass matrix components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      m11(iDelamEl)=m11(iDelamEl)+rho(i)*(h2(i)-h1(i));
      J11(iDelamEl)=J11(iDelamEl)+rho(i)*(h2(i)^3-h1(i)^3)/3;
    % stiffness matrix components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      a44(iDelamEl)=a44(iDelamEl)+s44(i)*5/4*(h2(i)-h1(i)-4/3*(h2(i)^3-h1(i)^3)/h^2);
      a45(iDelamEl)=a45(iDelamEl)+s45(i)*5/4*(h2(i)-h1(i)-4/3*(h2(i)^3-h1(i)^3)/h^2);
      a55(iDelamEl)=a55(iDelamEl)+s55(i)*5/4*(h2(i)-h1(i)-4/3*(h2(i)^3-h1(i)^3)/h^2);

      a11(iDelamEl)=a11(iDelamEl)+s11(i)*(h2(i)-h1(i));
      a12(iDelamEl)=a12(iDelamEl)+s12(i)*(h2(i)-h1(i));
      a16(iDelamEl)=a16(iDelamEl)+s16(i)*(h2(i)-h1(i));
      a22(iDelamEl)=a22(iDelamEl)+s22(i)*(h2(i)-h1(i));
      a26(iDelamEl)=a26(iDelamEl)+s26(i)*(h2(i)-h1(i));
      a66(iDelamEl)=a66(iDelamEl)+s66(i)*(h2(i)-h1(i));

      b11(iDelamEl)=b11(iDelamEl)+s11(i)*(h2(i)^2-h1(i)^2)/2;
      b12(iDelamEl)=b12(iDelamEl)+s12(i)*(h2(i)^2-h1(i)^2)/2;
      b16(iDelamEl)=b16(iDelamEl)+s16(i)*(h2(i)^2-h1(i)^2)/2;
      b22(iDelamEl)=b22(iDelamEl)+s22(i)*(h2(i)^2-h1(i)^2)/2;
      b26(iDelamEl)=b26(iDelamEl)+s26(i)*(h2(i)^2-h1(i)^2)/2;
      b66(iDelamEl)=b66(iDelamEl)+s66(i)*(h2(i)^2-h1(i)^2)/2;

      d11(iDelamEl)=d11(iDelamEl)+s11(i)*(h2(i)^3-h1(i)^3)/3;
      d12(iDelamEl)=d12(iDelamEl)+s12(i)*(h2(i)^3-h1(i)^3)/3;
      d16(iDelamEl)=d16(iDelamEl)+s16(i)*(h2(i)^3-h1(i)^3)/3;
      d22(iDelamEl)=d22(iDelamEl)+s22(i)*(h2(i)^3-h1(i)^3)/3;
      d26(iDelamEl)=d26(iDelamEl)+s26(i)*(h2(i)^3-h1(i)^3)/3;
      d66(iDelamEl)=d66(iDelamEl)+s66(i)*(h2(i)^3-h1(i)^3)/3;
    end
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

%% added mass layer 
addedMassEl=reshape(addedMassEl,1,[]);
imassEl = zeros(length(addedMassEl)*NofElNodes,1);
c=0;
for ne=addedMassEl
    c=c+1;
    n1=(ne-1)*NofElNodes+1;
    n2=n1+NofElNodes-1;
    I=n1:n2;
    imassEl((c-1)*NofElNodes+1:c*NofElNodes) = I; % indices of added mass nodes
end
% mass matrix components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m11(imassEl)=m11(imassEl)+rho_added_mass*(added_mass_thickness);
J11(imassEl)=J11(imassEl)+rho_added_mass*((h/2+added_mass_thickness/2)^3-(h/2)^3)/3;
 % stiffness matrix components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h2_=h/2+added_mass_thickness/2;
h1_=h/2; 
a44(imassEl)=a44(imassEl)+Qm(4,4)*5/4*(h2_-h1_-4/3*(h2_^3-h1_^3)/(h+added_mass_thickness)^2);
a45(imassEl)=a45(imassEl)+Qm(4,5)*5/4*(h2_-h1_-4/3*(h2_^3-h1_^3)/(h+added_mass_thickness)^2);
a55(imassEl)=a55(imassEl)+Qm(5,5)*5/4*(h2_-h1_-4/3*(h2_^3-h1_^3)/(h+added_mass_thickness)^2);

a11(imassEl)=a11(imassEl)+Qm(1,1)*(h2_-h1_);
a12(imassEl)=a12(imassEl)+Qm(1,2)*(h2_-h1_);
a16(imassEl)=a16(imassEl)+Qm(1,6)*(h2_-h1_);
a22(imassEl)=a22(imassEl)+Qm(2,2)*(h2_-h1_);
a26(imassEl)=a26(imassEl)+Qm(2,6)*(h2_-h1_);
a66(imassEl)=a66(imassEl)+Qm(6,6)*(h2_-h1_);

b11(imassEl)=b11(imassEl)+Qm(1,1)*(h2_^2-h1_^2)/2;
b12(imassEl)=b12(imassEl)+Qm(1,2)*(h2_^2-h1_^2)/2;
b16(imassEl)=b16(imassEl)+Qm(1,6)*(h2_^2-h1_^2)/2;
b22(imassEl)=b22(imassEl)+Qm(2,2)*(h2_^2-h1_^2)/2;
b26(imassEl)=b26(imassEl)+Qm(2,6)*(h2_^2-h1_^2)/2;
b66(imassEl)=b66(imassEl)+Qm(6,6)*(h2_^2-h1_^2)/2;

d11(imassEl)=d11(imassEl)+Qm(1,1)*(h2_^3-h1_^3)/3;
d12(imassEl)=d12(imassEl)+Qm(1,2)*(h2_^3-h1_^3)/3;
d16(imassEl)=d16(imassEl)+Qm(1,6)*(h2_^3-h1_^3)/3;
d22(imassEl)=d22(imassEl)+Qm(2,2)*(h2_^3-h1_^3)/3;
d26(imassEl)=d26(imassEl)+Qm(2,6)*(h2_^3-h1_^3)/3;
d66(imassEl)=d66(imassEl)+Qm(6,6)*(h2_^3-h1_^3)/3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of composite.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%