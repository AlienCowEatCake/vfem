x0 = 0;
y0 = 0;
z0 = 0;

x1 = 0.1;
y1 = 0.1;
z1 = 0.1;

cl = 0.025;

Point(1) = {x0, y0, z0, cl};
Point(2) = {x1, y0, z0, cl};
Point(3) = {x0, y1, z0, cl};
Point(4) = {x1, y1, z0, cl};

Point(5) = {x0, y0, z1, cl};
Point(6) = {x1, y0, z1, cl};
Point(7) = {x0, y1, z1, cl};
Point(8) = {x1, y1, z1, cl};

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

Surface Loop(25) = {14, 16, 18, 20, 22, 24};
Volume(26) = {25};

Physical Volume(27) = {26};

Physical Surface(28) = {14};
Physical Surface(29) = {16};
Physical Surface(30) = {18};
Physical Surface(31) = {20};
Physical Surface(32) = {22};
Physical Surface(33) = {24};

//Transfinite Surface "*" Alternate;
//Transfinite Volume "*";
