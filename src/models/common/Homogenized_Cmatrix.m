
rhom = 1200; % given 1150 - 1250 kg/m^3
rhof =1800; % given
em= 3.43e9; % assumed
ef=230e9; % given
nim=0.26; % assumed
nif=0.2; % assumed
vol=0.552; % given


[rho,e11,e22,e33,ni12,ni13,ni21,ni23,ni31,ni32,g12,g13,g23] = homogenization(rhom,rhof,em,ef,nim,nif,vol);
% something is wrong here
% shell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sierakowski page 46 eq. 2.33
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q11=e11./(1-ni12.*ni21);
q22=e22./(1-ni12.*ni21);
q44=g23;
q55=g13;
q66=g12;
q12=ni12.*e22./(1-ni12.*ni21);

[q11 q12 q22 q44 q55 q66 ]

%  1.0e+11 *
%
%   1.2907    0.0251    0.1107    0.0435    0.0451    0.0451
% rho = 1.5312e+03
% 1.2908    0.0259    0.1140    0.0448    0.0451  0.0451

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solid 3D
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

[q11 q12 q13 q22 q23 q33 q44 q55 q66]
% note that this is notation as in Sierakowski! = Voigt notation