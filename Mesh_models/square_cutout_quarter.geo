// Gmsh project created on Mon Nov 18 17:00:26 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 10, 0, 1.0};
//+
Point(2) = {0, 5, 0, 1.0};
//+
Point(3) = {5, 10, 0, 1.0};
//+
Point(4) = {4, 6, 0, 1.0};
//+
Point(5) = {4, 5, 0, 1.0};
//+
Point(6) = {5, 6, 0, 1.0};
//+
Line(1) = {1, 3};
//+
Line(2) = {3, 6};
//+
Line(3) = {6, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 2};
//+
Line(6) = {2, 1};
//+
Curve Loop(1) = {6, 1, 2, 3, 4, 5};
//+
Plane Surface(1) = {1};
//+
Physical Curve("Left", 7) = {6};
//+
Physical Curve("Top", 8) = {1};
//+
Physical Curve("Right", 9) = {2};
//+
Physical Curve("Bottom", 10) = {5};
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
//+
Transfinite Curve {6, 1} = 11 Using Progression 1;
//+
Transfinite Curve {3, 4} = 6 Using Progression 1;
//+
Physical Point("Left", 11) = {1, 2};
//+
Physical Point("Top", 12) = {1, 3};
//+
Physical Point("Right", 13) = {3, 6};
//+
Physical Point("Bottom", 14) = {2};
//+
Transfinite Curve {4, 3} = 6 Using Progression 1;
//+
Transfinite Curve {6, 1, 2, 5} = 11 Using Progression 1;
