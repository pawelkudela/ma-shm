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
Curve Loop(1) = {4, 1, 2, 3};
//+
Curve Loop(2) = {6};
//+
Curve Loop(3) = {5};
//+
Curve Loop(4) = {6};
//+
Curve Loop(5) = {5};
//+
Plane Surface(1) = {5};
//+
Plane Surface(2) = {4};
//+
Plane Surface(3) = {1, 2, 3};
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