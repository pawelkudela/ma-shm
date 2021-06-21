SetFactory("OpenCASCADE");
L=0.500000; 
W=0.500000; 
// delamination semi-major and semi-minor axis
a=0.01; 
b=0.005; 
// delamination centre coordinates
x=0.4; 
y=0.25; 
// delamination rotation
alpha=0*Pi; 
r=0.005000; 
xpzt=0.250000; 
ypzt=0.250000; 
S1x=0.250000; 
S1y=0.250000; 
S2x=0.35;
S2y=0.25;
S3x=0.25;
S3y=0.35;
S4x=0.32;
S4y=0.32;
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
Curve Loop(1) = {4, 1, 2, 3};
//+
Curve Loop(2) = {6};
//+
Curve Loop(3) = {5};
//+
Curve Loop(4) = {6};
//+
Curve Loop(5) = {5};
// pzt
Plane Surface(1) = {5};
// delam
Plane Surface(2) = {4};
// plate
Plane Surface(3) = {1, 2, 3};
//+
Point(9) = {S1x, S1y, 0, 1.0};
//+
Point{9} In Surface {1};// mesh node at S1
//+
Point(10) = {S2x, S2y, 0, 1.0};
//+
Point{10} In Surface {3};// mesh node at S2
//+
Point(11) = {S3x, S3y, 0, 1.0};
//+
Point{11} In Surface {3};// mesh node at S3
//+
Point(12) = {S4x, S4y, 0, 1.0};
//+
Point{12} In Surface {3};// mesh node at S4
//+
Physical Surface("pzt") = {1};
//+
Physical Surface("delam") = {2};
//+
Physical Surface("host") = {3};

Mesh.Algorithm = 6; // Frontal
Mesh.CharacteristicLengthFactor = 0.068000;
Mesh.CharacteristicLengthMin = 0.001000;
Mesh.CharacteristicLengthMax = 0.200000;
Mesh.RecombineAll = 1; // Apply recombination algorithm
Mesh.SubdivisionAlgorithm = 1; // 1: all quadrangles
Mesh.Smoothing = 1;
