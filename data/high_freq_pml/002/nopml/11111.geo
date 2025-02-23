cl__0 = 0.00009;
cl__1 = 0.00025;
cl__2 = 0.0005;
cl__3 = 0.0005;
//Point(1) = {0, 0, 0, 0.1};
//Point(2) = {0, 0.5, 0, 0.1};
//Point(3) = {0.5, 0, 0, 0.1};
//Point(4) = {-0.5, 0, 0, 0.1};
//Point(5) = {0, -0.5, 0, 0.1};
//Point(6) = {0, 0, 0.5, 0.1};
//Point(7) = {0, 0, -0.5, 0.1};
hh = 0.025 / 4;
hhPML = 0.025 / 3;
Point(8) = {-hh, -hh, -hh, cl__2};
Point(9) = {-hh, hh, -hh, cl__2};
Point(10) = {hh, -hh, -hh, cl__2};
Point(11) = {hh, hh, -hh, cl__2};
Point(12) = {-hh, -hh, hh, cl__2};
Point(13) = {-hh, hh, hh, cl__2};
Point(14) = {hh, -hh, hh, cl__2};
Point(15) = {hh, hh, hh, cl__2};
Point(16) = {-hhPML, -hhPML, -hhPML, cl__3};
Point(17) = {-hhPML, hhPML, -hhPML, cl__3};
Point(18) = {hhPML, -hhPML, -hhPML, cl__3};
Point(19) = {hhPML, hhPML, -hhPML, cl__3};
Point(20) = {-hhPML, -hhPML, hhPML, cl__3};
Point(21) = {-hhPML, hhPML, hhPML, cl__3};
Point(22) = {hhPML, -hhPML, hhPML, cl__3};
Point(23) = {hhPML, hhPML, hhPML, cl__3};

// Петля
za1 = 0.025 / 7;
diff = 0.0;
dCircle = 0.025 / 20;
Point(24) = {0, 0, za1, cl__0};
Point(25) = {dCircle, 0, za1, cl__0};
Point(26) = {0, dCircle, za1, cl__0};
Point(27) = {-dCircle, 0, za1, cl__0};
Point(28) = {0, -dCircle, za1, cl__0};
Circle(13) = {25, 24, 26};
Circle(14) = {26, 24, 27};
Circle(15) = {27, 24, 28};
Circle(16) = {28, 24, 25};
Line Loop(81) = {13, 14, 15, 16};
Plane Surface(1) = {81};
Physical Line(1) = {13, 16, 15, 14};


// Квадрат вокруг петли
hhCircle = 0.025 / 12;

Point(29) = {-hhCircle, -hhCircle, za1, cl__1};
Point(30) = {-hhCircle, hhCircle, za1, cl__1};
Point(31) = {hhCircle, hhCircle, za1, cl__1};
Point(32) = {hhCircle, -hhCircle, za1, cl__1};
Line(69) = {29, 30};
Line(70) = {30, 31};
Line(71) = {31, 32};
Line(72) = {32, 29};
Line Loop(82) = {69, 70, 71, 72};
Plane Surface(2) = {82, 81};

Point(33) = {-hhCircle, -hhCircle, za1-hhCircle, cl__1};
Point(34) = {-hhCircle, hhCircle, za1-hhCircle, cl__1};
Point(35) = {hhCircle, hhCircle, za1-hhCircle, cl__1};
Point(36) = {hhCircle, -hhCircle, za1-hhCircle, cl__1};
Point(37) = {-hhCircle, -hhCircle, za1+hhCircle, cl__1};
Point(38) = {-hhCircle, hhCircle, za1+hhCircle, cl__1};
Point(39) = {hhCircle, hhCircle, za1+hhCircle, cl__1};
Point(40) = {hhCircle, -hhCircle, za1+hhCircle, cl__1};


Line(83) = {33, 34};
Line(84) = {34, 35};
Line(85) = {35, 36};
Line(86) = {36, 33};
Line(87) = {37, 38};
Line(88) = {38, 39};
Line(89) = {39, 40};
Line(90) = {40, 37};
Line(91) = {33, 29};
Line(92) = {29, 37};
Line(93) = {34, 30};
Line(94) = {30, 38};
Line(95) = {35, 31};
Line(96) = {31, 39};
Line(97) = {36, 32};
Line(98) = {32, 40};


Line Loop(99) = {83, 84, 85, 86};
Plane Surface(100) = {99};
Line Loop(101) = {93, 70, -95, -84};
Plane Surface(102) = {101};
Line Loop(103) = {97, 72, -91, -86};
Plane Surface(104) = {103};
Line Loop(105) = {83, 93, -69, -91};
Plane Surface(106) = {105};
Line Loop(107) = {71, -97, -85, 95};
Plane Surface(108) = {107};
Line Loop(109) = {69, 94, -87, -92};
Plane Surface(110) = {109};
Line Loop(111) = {70, 96, -88, -94};
Plane Surface(112) = {111};
Line Loop(113) = {96, 89, -98, -71};
Plane Surface(114) = {113};
Line Loop(115) = {72, 92, -90, -98};
Plane Surface(116) = {115};
Line Loop(117) = {87, 88, 89, 90};
Plane Surface(118) = {117};

Surface Loop(119) = {100, 106, 102, 108, 104, 2, 1};
Volume(120) = {119};
Surface Loop(121) = {1, 2, 116, 110, 112, 114, 118};
Volume(122) = {121};



//Circle(1) = {5, 1, 3};
//Circle(2) = {2, 1, 3};
//Circle(3) = {2, 1, 4};
//Circle(4) = {4, 1, 5};
//Circle(5) = {6, 1, 5};
//Circle(6) = {5, 1, 7};
//Circle(7) = {6, 1, 2};
//Circle(8) = {2, 1, 7};
//Circle(9) = {6, 1, 4};
//Circle(10) = {6, 1, 3};
//Circle(11) = {4, 1, 7};
//Circle(12) = {7, 1, 3};
Line(31) = {12, 13};
Line(32) = {13, 9};
Line(33) = {9, 8};
Line(34) = {8, 12};
Line(35) = {13, 15};
Line(36) = {15, 14};
Line(37) = {14, 12};
Line(38) = {14, 10};
Line(39) = {10, 8};
Line(40) = {10, 11};
Line(41) = {11, 9};
Line(42) = {15, 11};
Line(57) = {23, 21};
Line(58) = {21, 17};
Line(59) = {17, 19};
Line(60) = {19, 23};
Line(61) = {23, 22};
Line(62) = {22, 18};
Line(63) = {18, 19};
Line(64) = {18, 16};
Line(65) = {16, 17};
Line(66) = {22, 20};
Line(67) = {20, 16};
Line(68) = {20, 21};
//Line Loop(14) = {9, 4, -5};
//Ruled Surface(14) = {14};
//Line Loop(16) = {7, 3, -9};
//Ruled Surface(16) = {16};
//Line Loop(18) = {10, -2, -7};
//Ruled Surface(18) = {18};
//Line Loop(20) = {5, 1, -10};
//Ruled Surface(20) = {20};
//Line Loop(22) = {6, 12, -1};
//Ruled Surface(22) = {22};
//Line Loop(24) = {12, -2, 8};
//Ruled Surface(24) = {24};
//Line Loop(26) = {11, -8, 3};
//Ruled Surface(26) = {26};
//Line Loop(28) = {6, -11, 4};
//Ruled Surface(28) = {28};
Line Loop(44) = {37, -34, -39, -38};
Ruled Surface(44) = {44};
Line Loop(46) = {36, 37, 31, 35};
Ruled Surface(46) = {46};
Line Loop(48) = {35, 42, 41, -32};
Ruled Surface(48) = {48};
Line Loop(50) = {31, 32, 33, 34};
Ruled Surface(50) = {50};
Line Loop(52) = {36, 38, 40, -42};
Ruled Surface(52) = {52};
Line Loop(54) = {41, 33, -39, 40};
Ruled Surface(54) = {54};
Line Loop(70) = {68, 58, -65, -67};
Ruled Surface(70) = {70};
Line Loop(72) = {66, 68, -57, 61};
Ruled Surface(72) = {72};
Line Loop(74) = {60, 61, 62, 63};
Ruled Surface(74) = {74};
Line Loop(76) = {63, -59, -65, -64};
Ruled Surface(76) = {76};
Line Loop(78) = {59, 60, 57, 58};
Ruled Surface(78) = {78};
Line Loop(80) = {62, 64, -67, -66};
Ruled Surface(80) = {80};
//Surface Loop(30) = {14, 16, 18, 20, 22, 28, 26, 24};
//Volume(30) = {30};
// Surface Loop(56) = {46, 52, 44, 50, 48, 54, 14, 16, 18, 20, 22, 28, 26, 24, 100, 106, 102, 108, 104, 116, 110, 112, 114, 118};
Surface Loop(56) = {46, 52, 44, 50, 48, 54, 100, 106, 102, 108, 104, 116, 110, 112, 114, 118};
Volume(56) = {56};
Surface Loop(82) = {80, 74, 78, 76, 70, 72, 46, 52, 44, 50, 48, 54};
Volume(82) = {82};
//, 100, 106, 102, 108, 104, 116, 110, 112, 114, 118

Physical Volume(1) = {120, 122, 56};
Physical Volume(2) = {82};
Physical Surface(1) = {70, 72, 74, 76, 80, 78};
