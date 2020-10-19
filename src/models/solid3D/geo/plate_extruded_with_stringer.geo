// Gmsh project created on Mon Oct 05 09:07:45 2020
SetFactory("OpenCASCADE");
div=30; // mesh divisions
thick_div=1; // thickness division
L=0.5; W=0.5; // plate dimensions
xpzt=0.25;
ypzt=0.25;
r=0.005; // pzt radius
hpzt=0.0002; //pzt thickness
h1=0.002; // plate thickness
h2=0.0015; // omega stringer thickness
h3=0.0365; // omega stringer total height
s=0.116; // stringer total width
sfb=0.0225; // stringer bottom foot width
sh=0.030; // stringer top part width
// b=h2*Cos(32*Pi/180)-Tan(32*Pi/180)*(h2-Sin(32*Pi/180)*h2);
b = (h2/Cos(58)-h2)/Tan(58);
bs=(h3-h2)*Cos(58*(Pi/180))/Sin(58*(Pi/180)); // stringer angled part width
sft=sfb-b; // stringer top foot width
L1=0.3215; // stringer position
//+
Point(1) = {L, L1, 0, 1.0};
//+
Point(2) = {L, L1+s-sfb, 0, 1.0};
//+
Point(3) = {L, W, h1, 1.0};
//+
Point(4) = {L, W, 0, 1.0};
//+
Point(5) = {L, L1, h1, 1.0};
//+
Point(6) = {L, L1+sfb, h1, 1.0};
//+
Point(7) = {L, L1+sft, h1+h2, 1.0};
//+
Point(8) = {L, L1, h1+h2, 1.0};
//+
Point(9) = {L, L1+s, h1, 1.0};
//+
Point(10) = {L, L1+s, h1+h2, 1.0};
//+
Point(11) = {L, L1+s-sft, h1+h2, 1.0};
//+
Point(12) = {L, L1+s-sfb, h1, 1.0};
//+
Point(13) = {L, L1+sft+bs, h1+h3, 1.0};
//+
Point(14) = {L, L1+s-sft-bs, h1+h3, 1.0};
//+
Point(15) = {L, L1+sfb+bs, h1+h3-h2, 1.0};
//+
Point(16) = {L, L1+s-sfb-bs, h1+h3-h2, 1.0};
//+
Point(17) = {L, L1+sfb, 0, 1.0};
//+
Point(18) = {L, L1+s, 0, 1.0};
//+
Line(1) = {1,5};
//+
Line(2) = {2,12};
//+
Line(3) = {3,9};
//+
Line(4) = {9,12};
//+
Line(5) = {12,6};
//+
Line(6) = {6,17};
//+
Line(7) = {4,3};
//+
Line(8) = {1,17};
//+
Line(9) = {7,8};
//+
Line(10) = {8,5};
//+
Line(11) = {9,10};
//+
Line(12) = {10,11};
//+
Line(13) = {13,7};
//+
Line(14) = {11,14};
//+
Line(15) = {14,13};
//+
Line(16) = {12,16};
//+
Line(17) = {16,15};
//+
Line(18) = {15,6};
//+
Line(19) = {6,7};
//+
Line(20) = {11,12};
//+
Line(21) = {13,15};
//+
Line(22) = {14,16};
//+
Line(23) = {5,6};
//+
Line(24) = {18,9};
//+
Line(25) = {17,2};
//+
Line(26) = {2,18};
//+
Line(27) = {18,4};
//+
//+
Curve Loop(1) = {1, 23, 6, -8};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {6, 25, 2, 5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {2, -4, -24, -26};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {24, -3, -7, -27};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {10, 23, 19, 9};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {19, -13, 21, 18};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {21, -17, -22, 15};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {22, -16, -20, 14};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {20, -4, 11, 12};
//+
Plane Surface(9) = {9};
//+
Extrude {-L, 0, 0} {
  Surface{1}; Surface{2}; Surface{3}; Surface{4}; Surface{5}; Surface{6}; Surface{7}; Surface{8}; Surface{9}; Layers{div}; Recombine;
}
//+
Transfinite Curve {10, 1, 19, 6, 2, 20, 24, 11, 7, 49, 44, 71, 68, 39, 64, 59, 54, 34, 30, 52, 21,22} = thick_div Using Progression 1;
//+
//Transfinite Curve {28, 29, 51, 53, 33, 31, 56, 58, 61, 63, 36, 38, 67, 43, 41, 70, 48, 46} = div Using Progression 1;
//+
Transfinite Curve {9, 23, 8, 13, 18, 15, 17, 14, 16, 12, 4, 55, 32, 35, 57, 60, 65, 62, 69, 66, 72, 42, 45, 26} = 4 Using Progression 1;
//+
Transfinite Curve {40, 37, 5, 25} = 8 Using Progression 1;
//+
Transfinite Curve {3, 27, 50, 47} = 7 Using Progression 1;
//+
Transfinite Surface {29} Alternated;
//+
Transfinite Surface {11} Alternated;
//+
Transfinite Surface {13} Alternated;
//+
Transfinite Surface {37} Alternated;
//+
Transfinite Surface {35} Alternated;
//+
Transfinite Surface {17} Alternated;
//+
Transfinite Surface {15} Alternated;
//+
Transfinite Surface {44} Alternated;
//+
Transfinite Surface {19} Alternated;
//+
Transfinite Surface {21} Alternated;
//+
Transfinite Surface {23} Alternated;
//+
Transfinite Surface {25} Alternated;
//+
Transfinite Volume{1} = {1, 17, 6, 5, 19, 22, 21, 20};
//+
Transfinite Volume{5} = {5, 6, 7, 8, 20, 21, 30, 29};
//+
Transfinite Volume{2} = {17, 2, 12, 6, 22, 23, 24, 21};
//+
Transfinite Volume{6} = {6, 15, 13, 7, 21, 32, 31, 30};
//+
Transfinite Volume{7} = {15, 16, 14, 13, 32, 33, 34, 31};
//+
Transfinite Volume{8} = {16, 12, 11, 14, 33, 24, 35, 34};
//+
Transfinite Volume{3} = {2, 18, 9, 12, 23, 26, 25, 24};
//+
Transfinite Volume{9} = {12, 9, 10, 11, 24, 25, 36, 35};
//+
Transfinite Volume{4} = {18, 4, 3, 9, 26, 28, 27, 25};
//+
p1=newp; Point(p1) = {0, 0, 0, 1.0};
p2=newp; Point(p2) = {L, 0, 0, 1.0}; 

l1 = newll; Line(l1) = {p1,p2};//+
Line(74) = {38, 1};
//+
Line(75) = {19, 37};
//+
l2 = newll;Curve Loop(l2) = {73, 74, 28, 75}; //plate
//+
c1 = newc; Circle(c1) = {xpzt, ypzt, 0, r, 0, 2*Pi};
//+
l3 = newll; Curve Loop(l3) = {c1}; //pzt//+
Transfinite Curve {77} = 13 Using Progression 1;
//+
Transfinite Curve {73,80} = div Using Progression 1;
//+
Transfinite Curve {75, 74} = 33 Using Progression 1;
//+
s1 = news; Plane Surface(s1) = {l2,l3}; //plate
s2 = news; Plane Surface(s2) = {l3}; //pzt bottom of the plate
//+
p3 = newp; Point(p3) = {xpzt, ypzt, 0, 1.0};
//+
Point{p3} In Surface {s2};// mesh node at centre of pzt
//+
Extrude {0, 0, h1} {
  Surface{s1}; Surface{s2}; Layers{1}; Recombine;
}
//+
Extrude {0, 0, hpzt} {
  Surface{87}; Layers{1}; Recombine;
}
Transfinite Curve {79, 78} = thick_div Using Progression 1;
//+

//+
Physical Volume("Plate1") = {10, 11};
//+
Physical Volume("Plate2") = {1};
//+
Physical Volume("Plate3") = {2};
//+
Physical Volume("Plate4") = {3};
//+
Physical Volume("Plate5") = {4};
//+
Physical Volume("Stringer_foot1") = {5};
//+
Physical Volume("Stringer_slope1") = {6};
//+
Physical Volume("Stringer_top") = {7};
//+
Physical Volume("Stringer_slope2") = {8};
//+
Physical Volume("Stringer_foot2") = {9};
//+
Physical Volume("pzt") = {12};
//+
// Mesh options
Mesh.RecombineAll = 1; // Apply recombination algorithm
Mesh.SubdivisionAlgorithm = 2; // 1: all quadrangles; 2: all hexas
Mesh.Smoothing = 1;
Mesh.CharacteristicLengthFactor = 0.08000;
Mesh.CharacteristicLengthMin = 0.001000;
Mesh.CharacteristicLengthMax = 0.200000;

Mesh.ElementOrder = 1;//+
// Options
General.SaveOptions = 0;
