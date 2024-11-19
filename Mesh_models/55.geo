// Gmsh project created on Mon Nov 18 15:33:05 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {5, 0, 0, 1.0};
//+
Point(2) = {5, 5, 0, 1.0};
//+
Point(3) = {0, 5, 0, 1.0};
//+
Point(4) = {0, 0, 0, 1.0};
//+
Line(1) = {4, 3};
//+
Line(2) = {3, 2};
//+
Line(3) = {2, 1};
//+
Line(4) = {1, 4};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Left", 5) = {1};
//+
Physical Curve("Top", 6) = {2};
//+
Physical Curve("Right", 7) = {3};
//+
Physical Curve("Bottom", 8) = {4};
//+
Physical Point("Bottom", 9) = {4, 1};
//+
Physical Point("Left", 10) = {4, 3};
//+
Physical Point("Top", 11) = {3, 2};
//+
Physical Point("Right", 12) = {2, 1};
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
//+
Transfinite Curve {1, 2, 3, 4} = 17 Using Progression 1;
//+
Transfinite Curve {1, 2, 3, 4} = 33 Using Progression 1;
//+
Transfinite Curve {1, 2, 3, 4} = 65 Using Progression 1;
//+
Transfinite Curve {1, 2, 3, 4} = 17 Using Progression 1;
//+
Transfinite Curve {1, 2, 3, 4} = 33 Using Progression 1;
//+
Transfinite Curve {1, 2, 3, 4} = 9 Using Progression 1;
