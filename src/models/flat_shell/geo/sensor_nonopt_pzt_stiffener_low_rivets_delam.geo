SetFactory("OpenCASCADE");
L=0.500000; 
W=0.500000; 
a=0.015000; 
b=0.01000; 
x=0.4000; 
y=0.100; 
alpha=0.000000 *Pi; 
r=0.005000; 
R=0.0025; // rivet radius
// PZT
xpzt1=0.185; 
ypzt1=0.45; 
xpzt2=0.185;
ypzt2=0.4;
xpzt3=0.185;
ypzt3=0.35;
xpzt4=0.185;
ypzt4=0.3;
xpzt5=0.185;
ypzt5=0.25;
xpzt6=0.185;
ypzt6=0.2;
xpzt7=0.185;
ypzt7=0.15;
xpzt8=0.185;
ypzt8=0.1;
xpzt9=0.185;
ypzt9=0.05;
xpzt10=0.315;
ypzt10=0.450;
xpzt11=0.315;
ypzt11=0.25;
xpzt12=0.315;
ypzt12=0.05;
xpzt13=0.3825;
ypzt13=0.450;
xpzt14=0.3825;
ypzt14=0.25;
xpzt15=0.3825;
ypzt15=0.05;
xpzt16=0.450;
ypzt16=0.450;
xpzt17=0.450;
ypzt17=0.250;
xpzt18=0.450;
ypzt18=0.05;
// Stiffener
xs1=0.235;
ys1=0.5;
xs2=0.235;
ys2=0;
xs3=0.265;
ys3=0.5;
xs4=0.265;
ys4=0;
// rivets
xr=0.25;
yr1=0.45;
yr2=0.4;
yr3=0.35;
yr4=0.3;
yr5=0.25;
yr6=0.2;
yr7=0.15;
yr8=0.1;
yr9=0.05;
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {L, 0, 0, 1.0};
//+
Point(3) = {L, W, 0, 1.0};
//+
Point(4) = {0, W, 0, 1.0};
//+
Point(5) = {xs1, ys1, 0, 1.0};
//+
Point(6) = {xs2, ys2, 0, 1.0};
//+
Point(7) = {xs3, ys3, 0, 1.0};
//+
Point(8) = {xs4, ys4, 0, 1.0};
//+
Point(9) = {xr, yr1, 0, 1.0};
//+
Point(10) = {xr, yr2, 0, 1.0};
//+
Point(11) = {xr, yr3, 0, 1.0};
//+
Point(12) = {xr, yr4, 0, 1.0};
//+
Point(13) = {xr, yr5, 0, 1.0};
//+
Point(14) = {xr, yr6, 0, 1.0};
//+
Point(15) = {xr, yr7, 0, 1.0};
//+
Point(16) = {xr, yr8, 0, 1.0};
//+
Point(17) = {xr, yr9, 0, 1.0};
//+
Line(1) = {1, 6};
//+
Line(2) = {6, 5};
//+
Line(3) = {5, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {6, 8};
//+
Line(6) = {8, 7};
//+
Line(7) = {7, 5};
//+
Line(8) = {5, 6};
//+
Line(9) = {8, 2};
//+
Line(10) = {2, 3};
//+
Line(11) = {3, 7};
//+
Line(12) = {7, 8};
//+
Circle(13) = {xpzt1, ypzt1, -0, r, 0, 2*Pi};
//+
Circle(14) = {xpzt2, ypzt2, -0, r, 0, 2*Pi};
//+
Circle(15) = {xpzt3, ypzt3, -0, r, 0, 2*Pi};
//+
Circle(16) = {xpzt4, ypzt4, -0, r, 0, 2*Pi};
//+
Circle(17) = {xpzt5, ypzt5, -0, r, 0, 2*Pi};
//+
Circle(18) = {xpzt6, ypzt6, -0, r, 0, 2*Pi};
//+
Circle(19) = {xpzt7, ypzt7, -0, r, 0, 2*Pi};
//+
Circle(20) = {xpzt8, ypzt8, -0, r, 0, 2*Pi};
//+
Ellipse(21) = {x, y, -0, a, b, 0, 2*Pi};
//+
Rotate {{0, 0, 1}, {x, y, 0}, alpha} {
  Curve{21}; 
}
Circle(22) = {xr, yr1, -0, r, 0, 2*Pi}; // rivet head
//+ 
Circle(23) = {xr, yr1, -0, R, 0, 2*Pi}; // rivet pin
//+ 
Circle(24) = {xr, yr2, -0, r, 0, 2*Pi}; // rivet head
//+ 
Circle(25) = {xr, yr2, -0, R, 0, 2*Pi}; // rivet pin
//+ 
Circle(26) = {xr, yr3, -0, r, 0, 2*Pi}; // rivet head
//+ 
Circle(27) = {xr, yr3, -0, R, 0, 2*Pi}; // rivet pin
//+ 
Circle(28) = {xr, yr4, -0, r, 0, 2*Pi}; // rivet head
//+ 
Circle(29) = {xr, yr4, -0, R, 0, 2*Pi}; // rivet pin
//+ 
Circle(30) = {xr, yr5, -0, r, 0, 2*Pi}; // rivet head
//+ 
Circle(31) = {xr, yr5, -0, R, 0, 2*Pi}; // rivet pin
//+ 
Circle(32) = {xr, yr6, -0, r, 0, 2*Pi}; // rivet head
//+ 
Circle(33) = {xr, yr6, -0, R, 0, 2*Pi}; // rivet pin
//+ 
Circle(34) = {xr, yr7, -0, r, 0, 2*Pi}; // rivet head
//+ 
Circle(35) = {xr, yr7, -0, R, 0, 2*Pi}; // rivet pin
//+ 
Circle(36) = {xr, yr8, -0, r, 0, 2*Pi}; // rivet head
//+ 
Circle(37) = {xr, yr8, -0, R, 0, 2*Pi}; // rivet pin
//+ 
Circle(38) = {xr, yr9, -0, r, 0, 2*Pi}; // rivet head
//+ 
Circle(39) = {xr, yr9, -0, R, 0, 2*Pi}; // rivet pin
//+ 
Circle(40) = {xpzt9, ypzt9, -0, r, 0, 2*Pi};
//+
Circle(41) = {xpzt10, ypzt10, -0, r, 0, 2*Pi};
//+
Circle(42) = {xpzt11, ypzt11, -0, r, 0, 2*Pi};
//+
Circle(43) = {xpzt12, ypzt12, -0, r, 0, 2*Pi};
//+
Circle(44) = {xpzt13, ypzt13, -0, r, 0, 2*Pi};
//+
Circle(45) = {xpzt14, ypzt14, -0, r, 0, 2*Pi};
//+
Circle(46) = {xpzt15, ypzt15, -0, r, 0, 2*Pi};
//+
Circle(47) = {xpzt16, ypzt16, -0, r, 0, 2*Pi};
//+
Circle(48) = {xpzt17, ypzt17, -0, r, 0, 2*Pi};
//+
Circle(49) = {xpzt18, ypzt18, -0, r, 0, 2*Pi};
//+
Curve Loop(1) = {1, 5,9,10,11,7,3,4}; //+ Plate
//+ 
Curve Loop(2) = {5, 6, 7, 8}; //+ stiffener
//+ 
Curve Loop(3) = {21}; //+ Delamination
//+
Curve Loop(4) = {13}; // PZT1
//+
Curve Loop(5) = {14}; // PZT2
//+
Curve Loop(6) = {15}; // PZT3
//+
Curve Loop(7) = {16}; // PZT4
//+
Curve Loop(8) = {17}; // PZT5
//+
Curve Loop(9) = {18}; // PZT6
//+
Curve Loop(10) = {19}; // PZT7
//+
Curve Loop(11) = {20}; // PZT8
//+
Curve Loop(12) = {22}; // Rivet head
//+
Curve Loop(13) = {23}; // Rivet pin 1
//+
Curve Loop(14) = {24}; // Rivet head
//+
Curve Loop(15) = {25}; // Rivet pin 2
//+
Curve Loop(16) = {26}; // Rivet head
//+
Curve Loop(17) = {27}; // Rivet pin 3
//+
Curve Loop(18) = {28}; // Rivet head
//+
Curve Loop(19) = {29}; // Rivet pin 4
//+
Curve Loop(20) = {30}; // Rivet head
//+ 
Curve Loop(21) = {31}; // Rivet pin 5
//+
Curve Loop(22) = {32}; // Rivet head 
//+  
Curve Loop(23) = {33}; // Rivet pin 6
//+
Curve Loop(24) = {34}; // Rivet head
//+
Curve Loop(25) = {35}; // Rivet pin 7
//+
Curve Loop(26) = {36}; // Rivet head
//+
Curve Loop(27) = {37}; // Rivet pin 8
//+
Curve Loop(28) = {38}; // Rivet head
//+
Curve Loop(29) = {39}; // Rivet pin 9
//+
Curve Loop(30) = {9,10,11,12}; // Right plate
//+
Curve Loop(31) = {40}; // PZT9
//+
Curve Loop(32) = {41}; // PZT10
//+
Curve Loop(33) = {42}; // PZT11
//+
Curve Loop(34) = {43}; // PZT12
//+
Curve Loop(35) = {44}; // PZT13
//+
Curve Loop(36) = {45}; // PZT14
//+
Curve Loop(37) = {46}; // PZT15
//+
Curve Loop(38) = {47}; // PZT16
//+
Curve Loop(39) = {48}; // PZT17
//+
Curve Loop(40) = {49}; // PZT18
//+
Plane Surface(1) = {4}; //PZT1
//+
Plane Surface(2) = {5}; //PZT2
//+
Plane Surface(3) = {6}; //PZT3
//+
Plane Surface(4) = {7}; //PZT4
//+
Plane Surface(5) = {8}; //PZT5
//+
Plane Surface(6) = {9}; //PZT6
//+
Plane Surface(7) = {10}; //PZT7
//+
Plane Surface(8) = {11}; //PZT8
//+
Plane Surface(9) = {3}; //delam
//+
Plane Surface(10) = {13}; //rivet pin 1
//+
Plane Surface(11) = {12,13}; //rivet head 1
//+
Plane Surface(12) = {15}; //rivet pin 2
//+
Plane Surface(13) = {14,15}; //rivet head 2
//+
Plane Surface(14) = {17}; //rivet pin 3
//+
Plane Surface(15) = {16,17}; //rivet head 3
//+
Plane Surface(16) = {19}; //rivet pin 4
//+
Plane Surface(17) = {18,19}; //rivet head 4
//+
Plane Surface(18) = {21}; //rivet pin 5
//+
Plane Surface(19) = {20,21}; //rivet head 5
//+
Plane Surface(20) = {23}; //rivet pin 6
//+
Plane Surface(21) = {22,23}; //rivet head 6
//+
Plane Surface(22) = {25}; //rivet pin 7
//+
Plane Surface(23) = {24,25}; //rivet head 7
//+
Plane Surface(24) = {27}; //rivet pin 8
//+
Plane Surface(25) = {26,27}; //rivet head 8
//+
Plane Surface(26) = {29}; //rivet pin 9
//+
Plane Surface(27) = {28,29}; //rivet head 9
//+
Plane Surface(28) = {2,12,14,16,18,20,22,24,26,28}; //stiffener without rivet heads
//+
Plane Surface(29) = {1, 3, 4, 5,6,7,8,9,10,11,31,2};//plate left
Plane Surface(30) = {30, 3,32,33,34,35,36,37,38,39,40};//plate right
//+
Plane Surface(31) = {31}; //PZT9
//+
Plane Surface(32) = {32}; //PZT10
//+
Plane Surface(33) = {33}; //PZT11
//+
Plane Surface(34) = {34}; //PZT12
//+
Plane Surface(35) = {35}; //PZT13
//+
Plane Surface(36) = {36}; //PZT14
//+
Plane Surface(37) = {37}; //PZT15
//+
Plane Surface(38) = {38}; //PZT16
//+
Plane Surface(39) = {39}; //PZT17
//+
Plane Surface(40) = {40}; //PZT18
//
Point(55) = {xpzt1, ypzt1, 0, 1.0};
//+
Point{55} In Surface {1};// mesh node at centre of pzt1
//+
Point(56) = {xpzt2, ypzt2, 0, 1.0};
//+
Point{56} In Surface {2};// mesh node at centre of pzt2
//+
Point(57) = {xpzt3, ypzt3, 0, 1.0};
//+
Point{57} In Surface {3};// mesh node at centre of pzt3
//+
Point(58) = {xpzt4, ypzt4, 0, 1.0};
//+
Point{58} In Surface {4};// mesh node at centre of pzt4
//+
Point(59) = {xpzt5, ypzt5, 0, 1.0};
//+
Point{59} In Surface {5};// mesh node at centre of pzt5
//+
Point(60) = {xpzt6, ypzt6, 0, 1.0};
//+
Point{60} In Surface {6};// mesh node at centre of pzt6
//+
Point(61) = {xpzt7, ypzt7, 0, 1.0};
//+
Point{61} In Surface {7};// mesh node at centre of pzt7
//+
Point(62) = {xpzt8, ypzt8, 0, 1.0};
//+
Point{62} In Surface {8};// mesh node at centre of pzt8
//+
Point(63) = {xpzt9, ypzt9, 0, 1.0};
//+
Point{63} In Surface {31};// mesh node at centre of pzt9
//+
Point(64) = {xpzt10, ypzt10, 0, 1.0};
//+
Point{64} In Surface {32};// mesh node at centre of pzt10
//+
Point(65) = {xpzt11, ypzt11, 0, 1.0};
//+
Point{65} In Surface {33};// mesh node at centre of pzt11
//+
Point(66) = {xpzt12, ypzt12, 0, 1.0};
//+
Point{66} In Surface {34};// mesh node at centre of pzt12
//+
Point(67) = {xpzt13, ypzt13, 0, 1.0};
//+
Point{67} In Surface {35};// mesh node at centre of pzt13
//+
Point(68) = {xpzt14, ypzt14, 0, 1.0};
//+
Point{68} In Surface {36};// mesh node at centre of pzt14
//+
Point(69) = {xpzt15, ypzt15, 0, 1.0};
//+
Point{69} In Surface {37};// mesh node at centre of pzt15
//+
Point(70) = {xpzt16, ypzt16, 0, 1.0};
//+
Point{70} In Surface {38};// mesh node at centre of pzt16
//+
Point(71) = {xpzt17, ypzt17, 0, 1.0};
//+
Point{71} In Surface {39};// mesh node at centre of pzt17
//+
Point(72) = {xpzt18, ypzt18, 0, 1.0};
//+
Point{72} In Surface {40};// mesh node at centre of pzt18
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
Physical Surface("pzt9") = {31};
//+
Physical Surface("pzt10") = {32};
//+
Physical Surface("pzt11") = {33};
//+
Physical Surface("pzt12") = {34};
//+
Physical Surface("pzt13") = {35};
//+
Physical Surface("pzt14") = {36};
//+
Physical Surface("pzt15") = {37};
//+
Physical Surface("pzt16") = {38};
//+
Physical Surface("pzt17") = {39};
//+
Physical Surface("pzt18") = {40};
//+
Physical Surface("delam") = {9};
//+
Physical Surface("rivet_pin1") = {10};
//+
Physical Surface("rivet_head1") = {11};
//+
Physical Surface("rivet_pin2") = {12};
//+
Physical Surface("rivet_head2") = {13};
//+
Physical Surface("rivet_pin3") = {14};
//+
Physical Surface("rivet_head_3") = {15};
//+
Physical Surface("rivet_pin4") = {16};
//+
Physical Surface("rivet_head4") = {17};
//+
Physical Surface("rivet_pin5") = {18};
//+
Physical Surface("rivet_head5") = {19};
//+
Physical Surface("rivet_pin6") = {20};
//+
Physical Surface("rivet_head6") = {21};
//+
Physical Surface("rivet_pin7") = {22};
//+
Physical Surface("rivet_head_7") = {23};
//+
Physical Surface("rivet_pin8") = {24};
//+
Physical Surface("rivet_head8") = {25};
//+
Physical Surface("rivet_pin9") = {26};
//+
Physical Surface("rivet_head9") = {27};
//+
Physical Surface("stiffener") = {28};
//+
Physical Surface("plate_left") = {29};
//+
Physical Surface("plate_right") = {30};
//+

Mesh.Algorithm = 6; // Frontal
Mesh.CharacteristicLengthFactor = 0.050000;
Mesh.CharacteristicLengthMin = 0.0005000;
Mesh.CharacteristicLengthMax = 0.200000;
Mesh.RecombineAll = 1; // Apply recombination algorithm
Mesh.SubdivisionAlgorithm = 1; // 1: all quadrangles
Mesh.Smoothing = 1;//+
