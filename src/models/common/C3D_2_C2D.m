function [q11,q12,q22,q44,q55,q66]=C3D_2_C2D(C11,C12,C13,C22,C23,C33,C44,C55,C66)

% convert 3D tensor C to 2D tensor C taking into account that epsz = 0 (for shell element)

C=[C11,C12,C13,   0,    0,    0;
     C12,C22,C23,   0,    0,    0;
     C13,C23,C33,   0,    0,    0;
       0,    0,    0, C44,   0 ,   0;
       0,    0,    0,    0,  C55,  0;
       0,    0,    0,    0,    0,  C66];

invC=inv(C);
e11=1/invC(1,1);
e22=1/invC(2,2);
ni12=-invC(2,1)/invC(1,1);
ni21=-invC(1,2)/invC(2,2);


q11=e11/(1-ni12*ni21);
q22=e22/(1-ni12*ni21);
q12= ni12*e22/(1-ni12*ni21);
q44=C44;
q55=C55;
q66=C66;

[q11 q12 q22 q44 q55 q66 ]


end