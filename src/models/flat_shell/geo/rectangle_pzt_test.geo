SetFactory("OpenCASCADE");
L=0.500000; 
W=0.500000; 
a=0.015000; 
b=0.010000; 
x=0.444444; 
y=0.444444; 
alpha=0.166667 *Pi; 
xpzt=0.250000; 
ypzt=0.250000;
Lpzt=0.02;
Wpzt=0.01; 
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {L, 0, 0, 1.0};
//+
Point(3) = {L, W, 0, 1.0};
//+
Point(4) = {0, W, 0, 1.0};
//pzt
Point(5) = {xpzt-Lpzt/2, ypzt-Wpzt/2,0,1};
//+
Point(6) = {xpzt+Lpzt/2, ypzt-Wpzt/2,0,1};
//+
Point(7) = {xpzt+Lpzt/2, ypzt+Wpzt/2,0,1};
//+
Point(8) = {xpzt-Lpzt/2, ypzt+Wpzt/2,0,1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {5, 6};
//+
Ellipse(6) = {x, y, -0, a, b, 0, 2*Pi};
//+
Rotate {{0, 0, 1}, {x, y, 0}, alpha} {
  Curve{6}; 
}
Line(7) = {6, 7};
//+
Line(8) = {7, 8};
//+
Line(9) = {8, 5};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Curve Loop(2) = {6};
//+
Curve Loop(3) = {5, 7, 8, 9};
//+
Plane Surface(1) = {3};
//+
Plane Surface(2) = {2};
//+
Plane Surface(3) = {1, 2, 3};
//+
Point(10) = {xpzt, ypzt, 0, 1.0};
//+
Point{10} In Surface {1};// mesh node at centre of pzt
//+
Physical Surface("pzt") = {1};
//+
Physical Surface("delam") = {2};
//+
Physical Surface("host") = {3};

Mesh.Algorithm = 6; // Frontal
Mesh.CharacteristicLengthFactor = 0.080000;
Mesh.CharacteristicLengthMin = 0.001000;
Mesh.CharacteristicLengthMax = 0.200000;
Mesh.RecombineAll = 1; // Apply recombination algorithm
Mesh.SubdivisionAlgorithm = 1; // 1: all quadrangles
Mesh.Smoothing = 1;