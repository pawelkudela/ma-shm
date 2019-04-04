SetFactory("OpenCASCADE");
a=0.020000; 
b=0.010000; 
x=0.000000; 
y=0.200000; 
alpha=0.666667 *Pi; 
r=0.005000; 
xpzt=0.250000; 
ypzt=0.250000; 
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.5, 0, 0, 1.0};
//+
Point(3) = {0.5, 0.5, 0, 1.0};
//+
Point(4) = {0, 0.5, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Circle(5) = {xpzt, ypzt, -0, r, 0, 2*Pi};
//+
Ellipse(6) = {x, y, -0, a, b, 0, 2*Pi};
//+
Rotate {{0, 0, 1}, {x, y, 0}, alpha} {
  Curve{6}; 
}
//+
BooleanDifference{ Curve{6}; Delete; }{ Curve{4}; }
//+
Delete {
  Curve{6}; 
}
Delete {
  Curve{8}; 
}
//+
Delete {
  Curve{4}; 
}
//+
Line(9) = {8, 7};
//+
Line(10) = {4, 8};
//+
Line(11) = {7, 1};
//+
Curve Loop(1) = {5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {9, 7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, -7, 11, 1, 2, 3};
//+
Curve Loop(4) = {5};
//+
Plane Surface(3) = {3, 4};
//+
Point(9) = {xpzt, ypzt, 0, 1.0};
//+
Point{9} In Surface {1};// mesh node at centre of pzt
//+
Physical Surface("pzt") = {1};
//+
Physical Surface("delam") = {2};
//+
Physical Surface("host") = {3};