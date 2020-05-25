function [D,J0]=composite3D_transf(h,h1,h2,rhom,rhof,em20,ef,nim,nif,vol,alpha,lay,Temp,Moist,prop,nx,ny,nz,Qz)


switch prop 
    case 'smeared'
        D=zeros(6,6*nx*ny);
    case 'layered'
        D=zeros(6,6*nx*ny*lay);
end
for cm=1:nx*ny
    for i=1:lay
    z=(h1(i)+h2(i))/h;
    for j=1:nz
        I(j)=(j-1)*nx*ny+cm;
    end
    N=shape1D_v2(nz,Qz,z);
    T=Temp(I);
    Ti=N*T;
    M=Moist(I);
    Mi=N*M;
    % modify Young's modulus of the matrix according to environmental conditions
    em(i)=environment_degradation(em20(i),Ti,Mi);
    end
% Note - function provided only for the same epoxy in each layer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mechanical properties of composite material components
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gm=em./(1+nim)/2; gf=ef./(1+nif)/2;
rho=rhof.*vol+rhom.*(1-vol)
e11=ef.*vol+em.*(1-vol)
e22=((ef+em)+(ef-em).*vol)./((ef+em)-(ef-em).*vol).*em
e33=e22;
ni12=nif.*vol+nim.*(1-vol)
ni13=ni12;
ni21=e22./e11.*ni12
ni31=e33./e11.*ni13
ni23=nif.*vol+nim.*(1-vol) .*(1.+nim-ni12.*em./e11) ./(1.-nim.^2 +nim.*ni12.*em./e11)
ni32=e33./e22.*ni23
g12=((gf+gm)+(gf-gm).*vol)./((gf+gm)-(gf-gm).*vol).*gm
g23=e22./2./(1+ni23)
g13=g12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sierakowski page 46 eq. 2.33
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Delta=(1-ni12.*ni21-ni13.*ni31-ni12.*ni23.*ni31-ni13.*ni21.*ni32-ni23.*ni32);
q11=e11.*(1-ni23.*ni32)./Delta;
q22=e22.*(1-ni31.*ni13)./Delta;
q33=e33.*(1-ni12.*ni21)./Delta;
q44=g23;
q55=g13;
q66=g12;
q12=(ni21+ ni31.*ni23).*e11./Delta;
q13=(ni31+ ni21.*ni32).*e11./Delta;
q23=(ni32+ ni12.*ni31).*e22./Delta;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% global single layer properties 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alphar=alpha*pi/180; m=cos(alphar); n=sin(alphar);
for k=1:lay
    if(abs(m(k))<eps) m(k)=0;end
    if(abs(n(k))<eps) n(k)=0;end
end
% elastic properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s11=q11.*m.^4+2*(q12+2*q66).*m.^2.*n.^2+q22.*n.^4;
s12=(q11+q22-4*q66).*m.^2.*n.^2+q12.*(m.^4+n.^4);
s13=q13.*m.^2 +q23.*n.^2;
s16=(q11-q12-2*q66).*m.^3.*n+(q12-q22+2*q66).*m.*n.^3;
s22=q11.*n.^4+2*(q12+2*q66).*m.^2.*n.^2+q22.*m.^4;
s23=n.^2.*q13 + m.^2.*q23;
s33=q33;
s26=(q11-q12-2*q66).*n.^3.*m+(q12-q22+2*q66).*n.*m.^3;
s36=(q13-q23).*m.*n;
s44=q44.*m.^2 + q55.*n.^2;
s45=(q55-q44).*m.*n;
s55=q55.*m.^2 + q44.*n.^2;
s66=(q11+q22-2*q12-2*q66).*m.^2.*n.^2+q66.*(m.^4+n.^4);

switch prop
    case 'smeared'
        J0=0; 
        a11=0; a12=0; a13=0; a16=0; a22=0; a23=0; a33=0; a26=0; a36=0; a44=0; a45=0; a55=0; a66=0; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % global multilayer properties
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i=1:lay;
            % mass matrix components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            J0=J0+rho(i)*(h2(i)-h1(i))/h;
            % stiffness matrix components %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            a11=a11+s11(i)*(h2(i)-h1(i))/h;
            a12=a12+s12(i)*(h2(i)-h1(i))/h;
            a13=a13+s13(i)*(h2(i)-h1(i))/h;
            a16=a16+s16(i)*(h2(i)-h1(i))/h;
            a22=a22+s22(i)*(h2(i)-h1(i))/h;
            a23=a23+s23(i)*(h2(i)-h1(i))/h;
            a33=a33+s33(i)*(h2(i)-h1(i))/h;
            a26=a26+s26(i)*(h2(i)-h1(i))/h;
            a36=a36+s36(i)*(h2(i)-h1(i))/h;
            a44=a44+s44(i)*(h2(i)-h1(i))/h;
            a45=a45+s45(i)*(h2(i)-h1(i))/h;
            a55=a55+s55(i)*(h2(i)-h1(i))/h;
            a66=a66+s66(i)*(h2(i)-h1(i))/h;
        end;
        D(1,cm*6-5)=a11;D(1,cm*6-4)=a12;D(1,cm*6-3)=a13;D(1,cm*6-2)=a16;
        D(2,cm*6-5)=a12;D(2,cm*6-4)=a22;D(2,cm*6-3)=a23;D(2,cm*6-2)=a26;
        D(3,cm*6-5)=a13;D(3,cm*6-4)=a23;D(3,cm*6-3)=a33;D(3,cm*6-2)=a36;
        D(4,cm*6-5)=a16;D(4,cm*6-4)=a26;D(4,cm*6-3)=a36;D(4,cm*6-2)=a66;
        D(5,cm*6-1)=a55;D(5,cm*6)=a45;
        D(6,cm*6-1)=a45;D(6,cm*6)=a44;
    case 'layered'
        J0=rho;
        for n=1:lay
            o=(n-1)*nx*ny*6;
        D(1,o+cm*6-5)=s11(n);D(1,o+cm*6-4)=s12(n);D(1,o+cm*6-3)=s13(n);D(1,o+cm*6-2)=s16(n);
        D(2,o+cm*6-5)=s12(n);D(2,o+cm*6-4)=s22(n);D(2,o+cm*6-3)=s23(n);D(2,o+cm*6-2)=s26(n);
        D(3,o+cm*6-5)=s13(n);D(3,o+cm*6-4)=s23(n);D(3,o+cm*6-3)=s33(n);D(3,o+cm*6-2)=s36(n);
        D(4,o+cm*6-5)=s16(n);D(4,o+cm*6-4)=s26(n);D(4,o+cm*6-3)=s36(n);D(4,o+cm*6-2)=s66(n);
        D(5,o+cm*6-1)=s55(n);D(5,o+cm*6)=s45(n);
        D(6,o+cm*6-1)=s45(n);D(6,o+cm*6)=s44(n);  
        end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of composite.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%