x1 = 10;
x2 = 50;
x3 = 700;

c1 = 2;
c2 = 10;
c3 = 100;

Point(1) = {x3, x3, x3, c3};
Point(2) = {x3, x3, -x3, c3};
Point(3) = {x3, -x3, x3, c3};
Point(4) = {x3, -x3, -x3, c3};
Point(5) = {-x3, x3, x3, c3};
Point(6) = {-x3, x3, -x3, c3};
Point(7) = {-x3, -x3, x3, c3};
Point(8) = {-x3, -x3, -x3, c3};
Line(1) = {4, 3};
Line(2) = {3, 1};
Line(3) = {1, 2};
Line(4) = {2, 4};
Line(5) = {4, 8};
Line(6) = {8, 6};
Line(7) = {6, 2};
Line(8) = {6, 5};
Line(9) = {5, 1};
Line(10) = {5, 7};
Line(11) = {7, 3};
Line(12) = {7, 8};

Point(9) = {0, 0, 0, c1};
Point(10) = {x1, 0, 0, c1};
Point(11) = {0, x1, 0, c1};
Point(12) = {-x1, 0, 0, c1};
Point(13) = {0, -x1, 0, c1};

Circle(13) = {10, 9, 11};
Circle(14) = {11, 9, 12};
Circle(15) = {12, 9, 13};
Circle(16) = {13, 9, 10};
Line Loop(17) = {13, 14, 15, 16};
Plane Surface(18) = {17};
Physical Line(19) = {13, 16, 15, 14};

Point(20) = {x2, x2, x2, c2};
Point(22) = {x2, x2, -x2, c2};
Point(23) = {x2, -x2, x2, c2};
Point(24) = {x2, -x2, -x2, c2};
Point(25) = {-x2, x2, x2, c2};
Point(26) = {-x2, x2, -x2, c2};
Point(27) = {-x2, -x2, x2, c2};
Point(28) = {-x2, -x2, -x2, c2};
Point(29) = {x2, x2, 0, c2};
Point(30) = {x2, -x2, 0, c2};
Point(31) = {-x2, x2, 0, c2};
Point(32) = {-x2, -x2, 0, c2};
Line(20) = {22, 29};
Line(21) = {29, 20};
Line(22) = {20, 23};
Line(23) = {23, 30};
Line(24) = {30, 24};
Line(25) = {24, 22};
Line(26) = {26, 31};
Line(27) = {31, 25};
Line(28) = {25, 27};
Line(29) = {27, 32};
Line(30) = {32, 28};
Line(31) = {28, 26};
Line(32) = {26, 22};
Line(33) = {24, 28};
Line(34) = {32, 31};
Line(35) = {30, 32};
Line(36) = {30, 29};
Line(37) = {29, 31};
Line(38) = {20, 25};
Line(39) = {23, 27};
Line Loop(40) = {36, 37, -34, -35};
Plane Surface(41) = {17, 40};
Line Loop(42) = {32, -25, 33, 31};
Plane Surface(43) = {42};
Line Loop(44) = {22, 39, -28, -38};
Plane Surface(45) = {44};
Line Loop(46) = {29, 34, 27, 28};
Plane Surface(47) = {46};
Line Loop(48) = {30, 31, 26, -34};
Plane Surface(49) = {48};
Line Loop(50) = {29, -35, -23, 39};
Plane Surface(51) = {50};
Line Loop(52) = {35, 30, -33, -24};
Plane Surface(53) = {52};
Line Loop(54) = {23, 36, 21, 22};
Plane Surface(55) = {54};
Line Loop(56) = {20, -36, 24, 25};
Plane Surface(57) = {56};
Line Loop(58) = {32, 20, 37, -26};
Plane Surface(59) = {58};
Line Loop(60) = {27, -38, -21, 37};
Plane Surface(61) = {60};
Surface Loop(62) = {43, 59, 57, 53, 49, 41, 18};
Volume(63) = {62};
Surface Loop(64) = {45, 55, 51, 47, 61, 41, 18};
Volume(65) = {64};

Line Loop(66) = {4, 5, 6, 7};
Plane Surface(67) = {66};
Line Loop(68) = {1, -11, 12, -5};
Plane Surface(69) = {68};
Line Loop(70) = {2, -9, 10, 11};
Plane Surface(71) = {70};
Line Loop(72) = {9, 3, -7, 8};
Plane Surface(73) = {72};
Line Loop(74) = {4, 1, 2, 3};
Plane Surface(75) = {74};
Line Loop(76) = {6, 8, 10, 12};
Plane Surface(77) = {76};
Surface Loop(78) = {67, 75, 69, 71, 73, 77};
Surface Loop(79) = {59, 43, 57, 53, 49, 61, 47, 51, 55, 45};
Volume(80) = {78, 79};
Physical Volume(81) = {63, 80, 65};

Physical Surface(82) = {67, 73, 71, 69, 77, 75};

