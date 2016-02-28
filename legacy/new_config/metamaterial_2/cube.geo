x0 = 0;
y0 = 0;
z0 = 0;

x1 = 0.003;
y1 = 0.003;
z1 = 0.003;

x2 = 0.0005;
y2 = 0.001;
z2 = 0.001; 

x3 = 0.001;
y3 = 0.002;
z3 = 0.002; 

x4 = 0.002;
y4 = 0.001;
z4 = 0.001; 

x5 = 0.0025;
y5 = 0.002;
z5 = 0.002;

cl1 = 0.0001;
cl2 = 0.00008;

Point(1) = {x0, y0, z0, cl1};
Point(2) = {x1, y0, z0, cl1};
Point(3) = {x0, y1, z0, cl1};
Point(4) = {x1, y1, z0, cl1};

Point(5) = {x0, y0, z1, cl1};
Point(6) = {x1, y0, z1, cl1};
Point(7) = {x0, y1, z1, cl1};
Point(8) = {x1, y1, z1, cl1};

Point(9) = {x2, y2, z2, cl2};
Point(10) = {x3, y2, z2, cl2};
Point(11) = {x2, y3, z2, cl2};
Point(12) = {x3, y3, z2, cl2};

Point(13) = {x2, y2, z3, cl2};
Point(14) = {x3, y2, z3, cl2};
Point(15) = {x2, y3, z3, cl2};
Point(66) = {x3, y3, z3, cl2};

Point(16) = {x4, y4, z4, cl2};
Point(17) = {x5, y4, z4, cl2};
Point(18) = {x4, y5, z4, cl2};
Point(19) = {x5, y5, z4, cl2};

Point(20) = {x4, y4, z5, cl2};
Point(21) = {x5, y4, z5, cl2};
Point(22) = {x4, y5, z5, cl2};
Point(23) = {x5, y5, z5, cl2};

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

Physical Surface(28) = {14};
Physical Surface(29) = {16};
Physical Surface(30) = {18};
Physical Surface(31) = {20};
Physical Surface(32) = {22};
Physical Surface(33) = {24};


Line(34) = {66, 12};
Line(35) = {12, 10};
Line(36) = {14, 10};
Line(37) = {66, 14};
Line(38) = {14, 13};
Line(39) = {10, 9};
Line(40) = {12, 11};
Line(41) = {66, 15};
Line(42) = {15, 11};
Line(43) = {11, 9};
Line(44) = {15, 13};
Line(45) = {13, 9};
Line Loop(46) = {42, -40, -34, 41};
Plane Surface(47) = {46};
Line Loop(48) = {35, 39, -43, -40};
Plane Surface(49) = {48};
Line Loop(50) = {36, 39, -45, -38};
Plane Surface(51) = {50};
Line Loop(52) = {43, -45, -44, 42};
Plane Surface(53) = {52};
Line Loop(54) = {35, -36, -37, 34};
Plane Surface(55) = {54};
Line Loop(56) = {37, 38, -44, -41};
Plane Surface(57) = {56};

Line(63) = {17, 16};
Line(64) = {16, 18};
Line(65) = {18, 19};
Line(66) = {19, 17};
Line(67) = {22, 20};
Line(68) = {20, 21};
Line(69) = {21, 23};
Line(70) = {23, 22};
Line(71) = {21, 17};
Line(72) = {20, 16};
Line(73) = {23, 19};
Line(74) = {22, 18};
Line Loop(75) = {63, -72, 68, 71};
Plane Surface(76) = {75};
Line Loop(77) = {71, -66, -73, -69};
Plane Surface(78) = {77};
Line Loop(79) = {72, 64, -74, 67};
Plane Surface(80) = {79};
Line Loop(81) = {68, 69, 70, 67};
Plane Surface(82) = {81};
Line Loop(83) = {73, -65, -74, -70};
Plane Surface(84) = {83};
Line Loop(85) = {66, 63, 64, 65};
Plane Surface(86) = {85};
Surface Loop(87) = {76, 86, 78, 84, 80, 82};
Volume(88) = {87};

Surface Loop(58) = {49, 55, 51, 53, 57, 47};
Volume(59) = {58};
Surface Loop(60) = {16, 24, 20, 18, 14, 22};
Volume(61) = {58, 87, 60};
Physical Volume(62) = {59};
Physical Volume(27) = {61};
Physical Volume(28) = {88};
