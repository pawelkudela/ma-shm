function [Sxx,Syy,Szz,Sxy,Syz,Sxz]=stresses_spec_p(Exx,Eyy,Ezz,Exy,Eyz,Exz,D11,D12,D13,D14,D22,D23,D24,D33,D34,D44,D55,D56,D66)

% calculate stresses at nodes of spectral elements
% parallel mode

Sxx=Exx.*D11+Eyy.*D12+Ezz.*D13+Exy.*D14;
Syy=Exx.*D12+Eyy.*D22+Ezz.*D23+Exy.*D24;
Szz=Exx.*D13+Eyy.*D23+Ezz.*D33+Exy.*D34;
Sxy=Exx.*D14+Eyy.*D24+Ezz.*D34+Exy.*D44;
Syz=Eyz.*D55+Exz.*D56;
Sxz=Eyz.*D56+Exz.*D66;