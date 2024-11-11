//+
Point(1) = {0, 10, 0, 1.0};
//+
Point(2) = {10, 10, 0, 1.0};
//+
Point(3) = {10, 0, 0, 1.0};
//+
Point(4) = {0, 0, 0, 1.0};
//+
Point(5) = {4, 4, 0, 1.0};
//+
Point(6) = {4, 6, 0, 1.0};
//+
Point(7) = {6, 6, 0, 1.0};
//+
Point(8) = {6, 4, 0, 1.0};
//+
Line(1) = {7, 2};
//+
Line(2) = {6, 1};
//+
Line(3) = {5, 4};
//+
Line(4) = {8, 3};
//+
Line(5) = {1, 2};
//+
Line(6) = {2, 3};
//+
Line(7) = {3, 4};
//+
Line(8) = {4, 1};
//+
Line(9) = {6, 7};
//+
Line(10) = {7, 8};
//+
Line(11) = {8, 5};
//+
Line(12) = {5, 6};
//+
Curve Loop(1) = {8, -2, -12, 3};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, -1, -9, 2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {6, -4, -10, 1};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {4, 7, -3, -11};
//+
Plane Surface(4) = {4};
//+
Physical Curve("Left", 13) = {8};
//+
Physical Curve("Top", 14) = {5};
//+
Physical Curve("Right", 15) = {6};
//+
Physical Curve("Bottom", 16) = {7};
//+
Physical Point("Bottom", 17) = {4, 3};
//+
Physical Point("Left", 18) = {4, 1};
//+
Physical Point("Top", 19) = {1, 2};
//+
Physical Point("Right", 20) = {2, 3};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Recombine Surface {1};
//+
Recombine Surface {2};
//+
Recombine Surface {3};
//+
Recombine Surface {4};
//+
Transfinite Curve {8, 5, 6, 7} = 11 Using Progression 1;
//+
Transfinite Curve {2, 1, 4, 3, 12, 9, 10, 11} = 11 Using Progression 1;
//+
Transfinite Curve {3, 2, 1, 4} = 6 Using Progression 1;
