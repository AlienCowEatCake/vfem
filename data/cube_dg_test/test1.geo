x0 = 0;
y0 = 0;
z0 = 0;

x1 = 0.5;

x2 = 1.0;
y2 = 1.0;
z2 = 1.0;

cl = 0.25;

Point(1) = {x0, y0, z0, cl};
Point(2) = {x1, y0, z0, cl};
Point(3) = {x0, y2, z0, cl};
Point(4) = {x1, y2, z0, cl};

Point(5) = {x0, y0, z2, cl};
Point(6) = {x1, y0, z2, cl};
Point(7) = {x0, y2, z2, cl};
Point(8) = {x1, y2, z2, cl};

Line(1) = {1, 2};
Line(2) = {1, 3};
Line(3) = {1, 5};
Line(4) = {2, 4};
Line(5) = {2, 6};
Line(6) = {3, 4};
Line(7) = {3, 7};
Line(8) = {4, 8};
Line(9) = {5, 6};
Line(10) = {5, 7};
Line(11) = {6, 8};
Line(12) = {7, 8};

Line Loop(13) = {1, 4, -6, -2};
Plane Surface(14) = {13};
Line Loop(15) = {8, -11, -5, 4};
Plane Surface(16) = {15};
Line Loop(17) = {2, 7, -10, -3};
Plane Surface(18) = {17};
Line Loop(19) = {11, -12, -10, 9};
Plane Surface(20) = {19};
Line Loop(21) = {5, -9, -3, 1};
Plane Surface(22) = {21};
Line Loop(23) = {8, -12, -7, 6};
Plane Surface(24) = {23};

Transfinite Surface(16) Alternate;

Surface Loop(25) = {14, 16, 18, 20, 22, 24};
Volume(26) = {25};

Point(11) = {x2, y0, z0, cl};
Point(12) = {x1, y0, z0, cl};
Point(13) = {x2, y2, z0, cl};
Point(14) = {x1, y2, z0, cl};

Point(15) = {x2, y0, z2, cl};
Point(16) = {x1, y0, z2, cl};
Point(17) = {x2, y2, z2, cl};
Point(18) = {x1, y2, z2, cl};

Line(21) = {11, 12};
Line(22) = {11, 13};
Line(23) = {11, 15};
Line(24) = {12, 14};
Line(25) = {12, 16};
Line(26) = {13, 14};
Line(27) = {13, 17};
Line(28) = {14, 18};
Line(29) = {15, 16};
Line(30) = {15, 17};
Line(31) = {16, 18};
Line(32) = {17, 18};

Line Loop(33) = {25, 31, -28, -24};
Plane Surface(34) = {33};
Line Loop(35) = {29, 31, -32, -30};
Plane Surface(36) = {35};
Line Loop(37) = {25, -29, -23, 21};
Plane Surface(38) = {37};
Line Loop(39) = {24, -26, -22, 21};
Plane Surface(40) = {39};
Line Loop(41) = {28, -32, -27, 26};
Plane Surface(42) = {41};
Line Loop(43) = {27, -30, -23, 22};
Plane Surface(44) = {43};

Surface Loop(45) = {34, 38, 36, 42, 44, 40};
Volume(46) = {45};

Transfinite Surface(34) Alternate;

Physical Volume(27) = {26, 46};

Physical Surface(28) = {14, 40};
Physical Surface(29) = {44};
Physical Surface(30) = {18};
Physical Surface(31) = {20, 36};
Physical Surface(32) = {22, 38};
Physical Surface(33) = {24, 42};

Physical Surface(41) = {16};
Physical Surface(42) = {34};
