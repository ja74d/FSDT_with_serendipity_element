// Gmsh project created on Sun Nov 17 10:06:13 2024
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 10, 0, 1.0};
//+
Point(2) = {10, 10, 0, 1.0};
//+
Point(3) = {10, 0, 0, 1.0};
//+
Point(4) = {0, 0, 0, 1.0};
//+
Point(5) = {4, 4.5, 0, 1.0};
//+
Point(6) = {5, 4, 0, 1.0};
//+
Point(7) = {4, 5.5, 0, 1.0};
//+
Point(8) = {5, 6, 0, 1.0};
//+
Point(9) = {6, 4.5, 0, 1.0};
//+
Point(10) = {6, 5.5, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {7, 8};
//+
Line(6) = {8, 10};
//+
Line(7) = {10, 9};
//+
Line(8) = {9, 6};
//+
Line(9) = {6, 5};
//+
Line(10) = {5, 7};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Curve Loop(2) = {5, 6, 7, 8, 9, 10};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("Left", 11) = {4};
//+
Physical Curve("Top", 12) = {1};
//+
Physical Curve("Right", 13) = {2};
//+
Physical Curve("Bottom", 14) = {3};
//+
Physical Point("Bottom", 15) = {4, 3};
//+
Physical Point("Left", 16) = {4, 1};
//+
Physical Point("Top", 17) = {1, 2};
//+
Physical Point("Right", 18) = {2, 3};
//+
Transfinite Surface {1};
//+
Recombine Surface {1};
//+
Transfinite Curve {4, 1, 2, 3} = 17 Using Progression 1;
//+
Transfinite Curve {5, 10, 9, 8, 7, 6} = 5 Using Progression 1;
