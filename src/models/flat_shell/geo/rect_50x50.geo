SetFactory("OpenCASCADE");
L=0.500000; 
W=0.500000; 
Nx=51;
Ny=51;
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
Curve Loop(1) = {4, 1, 2, 3};
//+

Plane Surface(1) = {1};

Physical Surface("plate") = {1};
//+
Transfinite Surface {1};
//+
Transfinite Curve {4, 2} = Nx Using Progression 1;
//+
Transfinite Curve {3, 1} = Ny Using Progression 1;
//+
//Mesh.Algorithm = 6; // Frontal
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.CharacteristicLengthFactor = 0.0;
Mesh.CharacteristicLengthMin = 0.00000;
Mesh.CharacteristicLengthMax = 0.00000;
Recombine Surface {1};
//Mesh.RecombineAll = 1; // Apply recombination algorithm
//Mesh.SubdivisionAlgorithm = 1; // 1: all quadrangles
//Mesh.Smoothing = 1;


