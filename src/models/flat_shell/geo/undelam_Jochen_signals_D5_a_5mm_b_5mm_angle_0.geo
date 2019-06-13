SetFactory("OpenCASCADE");
L=0.500000; 
W=0.500000; 
a=0.005000; 
b=0.005000; 
x=0.25000; 
y=0.42700; 
alpha=0.000000 *Pi; 
r=0.005000; 
xpzt1=0.450000; 
ypzt1=0.470000; 
xpzt2=0.370;
ypzt2=0.470;
xpzt3=0.290;
ypzt3=0.470;
xpzt4=0.210;
ypzt4=0.470;
xpzt5=0.130;
ypzt5=0.470;
xpzt6=0.050;
ypzt6=0.470;
xpzt7=0.450;
ypzt7=0.030;
xpzt8=0.370;
ypzt8=0.030;
xpzt9=0.290;
ypzt9=0.030;
xpzt10=0.210;
ypzt10=0.030;
xpzt11=0.130;
ypzt11=0.030;
xpzt12=0.050;
ypzt12=0.030;
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
Circle(5) = {xpzt1, ypzt1, -0, r, 0, 2*Pi};
//+
Circle(6) = {xpzt2, ypzt2, -0, r, 0, 2*Pi};
//+
Circle(7) = {xpzt3, ypzt3, -0, r, 0, 2*Pi};
//+
Circle(8) = {xpzt4, ypzt4, -0, r, 0, 2*Pi};
//+
Circle(9) = {xpzt5, ypzt5, -0, r, 0, 2*Pi};
//+
Circle(10) = {xpzt6, ypzt6, -0, r, 0, 2*Pi};
//+
Circle(11) = {xpzt7, ypzt7, -0, r, 0, 2*Pi};
//+
Circle(12) = {xpzt8, ypzt8, -0, r, 0, 2*Pi};
//+
Circle(13) = {xpzt9, ypzt9, -0, r, 0, 2*Pi};
//+
Circle(14) = {xpzt10, ypzt10, -0, r, 0, 2*Pi};
//+
Circle(15) = {xpzt11, ypzt11, -0, r, 0, 2*Pi};
//+
Circle(16) = {xpzt12, ypzt12, -0, r, 0, 2*Pi};
//+
Ellipse(17) = {x, y, -0, a, b, 0, 2*Pi};
//+
Rotate {{0, 0, 1}, {x, y, 0}, alpha} {
  Curve{17}; 
}
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Curve Loop(2) = {17};
//+
Curve Loop(3) = {5};
//+
Curve Loop(4) = {17};
//+
Curve Loop(5) = {5};
//+
Curve Loop(6) = {6};
//+
Curve Loop(7) = {7};
//+
Curve Loop(8) = {8};
//+
Curve Loop(9) = {9};
//+
Curve Loop(10) = {10};
//+
Curve Loop(11) = {11};
//+
Curve Loop(12) = {12};
//+
Curve Loop(13) = {13};
//+
Curve Loop(14) = {14};
//+
Curve Loop(15) = {15};
//+
Curve Loop(16) = {16};
//+
Plane Surface(1) = {5}; //PZT1
//+
Plane Surface(2) = {6}; //PZT2
//+
Plane Surface(3) = {7}; //PZT3
//+
Plane Surface(4) = {8}; //PZT4
//+
Plane Surface(5) = {9}; //PZT5
//+
Plane Surface(6) = {10}; //PZT6
//+
Plane Surface(7) = {11}; //PZT7
//+
Plane Surface(8) = {12}; //PZT8
//+
Plane Surface(9) = {13}; //PZT9
//+
Plane Surface(10) = {14}; //PZT10
//+
Plane Surface(11) = {15}; //PZT11
//+
Plane Surface(12) = {16}; //PZT12
//+
Plane Surface(13) = {4}; //delam
//+
Plane Surface(14) = {1, 2, 5,6,7,8,9,10,11,12,13,14,15,16};//plate+delam+PZT
//+
Point(29) = {xpzt1, ypzt1, 0, 1.0};
//+
Point{29} In Surface {1};// mesh node at centre of pzt1
//+
Point(18) = {xpzt2, ypzt2, 0, 1.0};
//+
Point{18} In Surface {2};// mesh node at centre of pzt2
//+
Point(19) = {xpzt3, ypzt3, 0, 1.0};
//+
Point{19} In Surface {3};// mesh node at centre of pzt3
//+
Point(20) = {xpzt4, ypzt4, 0, 1.0};
//+
Point{20} In Surface {4};// mesh node at centre of pzt4
//+
Point(21) = {xpzt5, ypzt5, 0, 1.0};
//+
Point{21} In Surface {5};// mesh node at centre of pzt5
//+
Point(22) = {xpzt6, ypzt6, 0, 1.0};
//+
Point{22} In Surface {6};// mesh node at centre of pzt6
//+
Point(23) = {xpzt7, ypzt7, 0, 1.0};
//+
Point{23} In Surface {7};// mesh node at centre of pzt7
//+
Point(24) = {xpzt8, ypzt8, 0, 1.0};
//+
Point{24} In Surface {8};// mesh node at centre of pzt8
//+
Point(25) = {xpzt9, ypzt9, 0, 1.0};
//+
Point{25} In Surface {9};// mesh node at centre of pzt9
//+
Point(26) = {xpzt10, ypzt10, 0, 1.0};
//+
Point{26} In Surface {10};// mesh node at centre of pzt10
//+
Point(27) = {xpzt11, ypzt11, 0, 1.0};
//+
Point{27} In Surface {11};// mesh node at centre of pzt11
//+
Point(28) = {xpzt12, ypzt12, 0, 1.0};
//+
Point{28} In Surface {12};// mesh node at centre of pzt12
//+
Physical Surface("pzt1") = {1};
//+
Physical Surface("pzt2") = {2};
//+
Physical Surface("pzt3") = {3};
//+
Physical Surface("pzt4") = {4};
//+
Physical Surface("pzt5") = {5};
//+
Physical Surface("pzt6") = {6};
//+
Physical Surface("pzt7") = {7};
//+
Physical Surface("pzt8") = {8};
//+
Physical Surface("pzt9") = {9};
//+
Physical Surface("pzt10") = {10};
//+
Physical Surface("pzt11") = {11};
//+
Physical Surface("pzt12") = {12};
//+
Physical Surface("delam") = {13};
//+
Physical Surface("host") = {14};

Mesh.Algorithm = 6; // Frontal
Mesh.CharacteristicLengthFactor = 0.080000;
Mesh.CharacteristicLengthMin = 0.001000;
Mesh.CharacteristicLengthMax = 0.200000;
Mesh.RecombineAll = 1; // Apply recombination algorithm
Mesh.SubdivisionAlgorithm = 1; // 1: all quadrangles
Mesh.Smoothing = 1;