//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {L, 0, 0, 1.0};
//+
Point(3) = {L, W, 0, 1.0};
//+
Point(4) = {0, W, 0, 1.0};
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
BooleanDifference{ Curve{6}; Delete; }{ Curve{1};Curve{4};}
//+
Delete {
  Curve{7}; 
}
//+
Delete {
  Curve{1}; 
}
//+
Delete {
  Curve{4}; 
}
//+
Line(9) = {1, 6};
//+
Line(10) = {7, 1};
//+
Line(11) = {6, 2};
//+
Line(12) = {4, 7};
//+
Curve Loop(1) = {5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {10,9, 6};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {6, 11, 2, 3, 12};
 
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