xArea = 2000;
xPmlEnd = 700;
xPmlBeg = 600;
xSource = 10;
x0 = 0;

zAirArea = 2000;
zAirPmlEnd = 700;
zAirPmlBeg = 600;
zSource = 10;
zWater = 0;
zWaterPmlBeg = -600;
zWaterPmlEnd = -700;
zWaterArea = -800;

kSource = 1.0;
kPmlBeg = 90.0;
kPmlEnd = 110.0;
kPmlBegC = 40.0;
kPmlEndC = 60.0;

// Источник
Point(1) = {xSource, xSource, zSource, kSource};
Point(2) = {xSource, -xSource, zSource, kSource};
Point(3) = {-xSource, xSource, zSource, kSource};
Point(4) = {-xSource, -xSource, zSource, kSource};
Point(5) = {xSource, xSource, zWater, kSource};
Point(6) = {xSource, -xSource, zWater, kSource};
Point(7) = {-xSource, xSource, zWater, kSource};
Point(8) = {-xSource, -xSource, zWater, kSource};
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 3};
Line(4) = {3, 1};
Line(5) = {3, 7};
Line(6) = {7, 8};
Line(7) = {8, 4};
Line(8) = {8, 6};
Line(9) = {6, 2};
Line(10) = {1, 5};
Line(11) = {5, 6};
Line(12) = {5, 7};
Line Loop(13) = {2, 3, 4, 1};
Plane Surface(14) = {13};
Line Loop(15) = {10, 11, 9, -1};
Plane Surface(16) = {15};
Line Loop(17) = {12, -5, 4, 10};
Plane Surface(18) = {17};
Line Loop(19) = {12, 6, 8, -11};
Plane Surface(20) = {19};
Line Loop(21) = {9, 2, -7, 8};
Plane Surface(22) = {21};
Line Loop(23) = {6, 7, 3, 5};
Plane Surface(24) = {23};
Surface Loop(25) = {14, 22, 16, 18, 20, 24};
Volume(26) = {25};

// Поверхность воды до PML
Point(9) = {xPmlBeg, xPmlBeg, zWater, kPmlBegC};
Point(10) = {xPmlBeg, -xPmlBeg, zWater, kPmlBegC};
Point(11) = {-xPmlBeg, xPmlBeg, zWater, kPmlBegC};
Point(12) = {-xPmlBeg, -xPmlBeg, zWater, kPmlBegC};
Line(27) = {11, 9};
Line(28) = {9, 10};
Line(29) = {10, 12};
Line(30) = {12, 11};
Line Loop(31) = {27, 28, 29, 30};
Plane Surface(32) = {19, 31};

// Воздух до PML
Point(13) = {xPmlBeg, xPmlBeg, zAirPmlBeg, kPmlBeg};
Point(14) = {xPmlBeg, -xPmlBeg, zAirPmlBeg, kPmlBeg};
Point(15) = {-xPmlBeg, xPmlBeg, zAirPmlBeg, kPmlBeg};
Point(16) = {-xPmlBeg, -xPmlBeg, zAirPmlBeg, kPmlBeg};
Line(33) = {11, 15};
Line(34) = {15, 13};
Line(35) = {13, 9};
Line(36) = {13, 14};
Line(37) = {14, 10};
Line(38) = {12, 16};
Line(39) = {16, 14};
Line(40) = {16, 15};
Line Loop(41) = {28, -37, -36, 35};
Plane Surface(42) = {41};
Line Loop(43) = {29, 38, 39, 37};
Plane Surface(44) = {43};
Line Loop(45) = {27, -35, -34, -33};
Plane Surface(46) = {45};
Line Loop(47) = {30, 33, -40, -38};
Plane Surface(48) = {47};
Line Loop(49) = {39, -36, -34, -40};
Plane Surface(50) = {49};
Surface Loop(51) = {50, 44, 32, 48, 46, 42, 14, 22, 16, 18, 24};
Volume(52) = {51};

// Поверхность воды в PML
Point(17) = {xPmlEnd, xPmlEnd, zWater, kPmlEndC};
Point(18) = {xPmlEnd, -xPmlEnd, zWater, kPmlEndC};
Point(19) = {-xPmlEnd, xPmlEnd, zWater, kPmlEndC};
Point(20) = {-xPmlEnd, -xPmlEnd, zWater, kPmlEndC};
Line(53) = {19, 20};
Line(54) = {20, 18};
Line(55) = {18, 17};
Line(56) = {17, 19};
Line Loop(57) = {56, 53, 54, 55};
Plane Surface(58) = {31, 57};

// Воздух в PML
Point(21) = {xPmlEnd, xPmlEnd, zAirPmlEnd, kPmlBeg};
Point(22) = {xPmlEnd, -xPmlEnd, zAirPmlEnd, kPmlBeg};
Point(23) = {-xPmlEnd, xPmlEnd, zAirPmlEnd, kPmlBeg};
Point(24) = {-xPmlEnd, -xPmlEnd, zAirPmlEnd, kPmlBeg};
Line(59) = {19, 23};
Line(60) = {20, 24};
Line(61) = {18, 22};
Line(62) = {17, 21};
Line(63) = {21, 22};
Line(64) = {21, 23};
Line(65) = {23, 24};
Line(66) = {24, 22};
Line Loop(67) = {55, 62, 63, -61};
Plane Surface(68) = {67};
Line Loop(69) = {54, 61, -66, -60};
Plane Surface(70) = {69};
Line Loop(71) = {53, 60, -65, -59};
Plane Surface(72) = {71};
Line Loop(73) = {56, 59, -64, -62};
Plane Surface(74) = {73};
Line Loop(75) = {66, -63, 64, 65};
Plane Surface(76) = {75};
Surface Loop(77) = {76, 70, 58, 68, 74, 72, 48, 46, 42, 44, 50};
Volume(78) = {77};

// Вода до PML
Point(25) = {xPmlBeg, xPmlBeg, zWaterPmlBeg, kPmlBeg};
Point(26) = {xPmlBeg, -xPmlBeg, zWaterPmlBeg, kPmlBeg};
Point(27) = {-xPmlBeg, xPmlBeg, zWaterPmlBeg, kPmlBeg};
Point(28) = {-xPmlBeg, -xPmlBeg, zWaterPmlBeg, kPmlBeg};
Line(79) = {28, 26};
Line(80) = {26, 25};
Line(81) = {25, 27};
Line(82) = {27, 28};
Line(83) = {10, 26};
Line(84) = {9, 25};
Line(85) = {11, 27};
Line(86) = {12, 28};
Line Loop(87) = {29, 86, 79, -83};
Plane Surface(88) = {87};
Line Loop(89) = {83, 80, -84, 28};
Plane Surface(90) = {89};
Line Loop(91) = {84, 81, -85, 27};
Plane Surface(92) = {91};
Line Loop(93) = {30, 85, 82, -86};
Plane Surface(94) = {93};
Line Loop(95) = {79, 80, 81, 82};
Plane Surface(96) = {95};
Surface Loop(97) = {96, 88, 94, 92, 90, 32, 20};
Volume(98) = {97};

// Вода в PML
Point(29) = {xPmlEnd, xPmlEnd, zWaterPmlEnd, kPmlEnd};
Point(30) = {xPmlEnd, -xPmlEnd, zWaterPmlEnd, kPmlEnd};
Point(31) = {-xPmlEnd, xPmlEnd, zWaterPmlEnd, kPmlEnd};
Point(32) = {-xPmlEnd, -xPmlEnd, zWaterPmlEnd, kPmlEnd};
Line(99) = {29, 31};
Line(100) = {31, 32};
Line(101) = {32, 30};
Line(102) = {30, 29};
Line(103) = {29, 17};
Line(104) = {30, 18};
Line(105) = {31, 19};
Line(106) = {32, 20};
Line Loop(107) = {100, 101, 102, 99};
Plane Surface(108) = {107};
Line Loop(109) = {56, -105, -99, 103};
Plane Surface(110) = {109};
Line Loop(111) = {103, -55, -104, 102};
Plane Surface(112) = {111};
Line Loop(113) = {101, 104, -54, -106};
Plane Surface(114) = {113};
Line Loop(115) = {100, 106, -53, -105};
Plane Surface(116) = {115};
Surface Loop(117) = {116, 108, 114, 112, 110, 58, 96, 88, 94, 92, 90};
Volume(118) = {117};

// Воздух
Physical Volume(119) = {52, 26};
// Воздух-PML
Physical Volume(120) = {78};
// Вода
Physical Volume(121) = {98};
// Вода-PML
Physical Volume(122) = {118};
// Первые краевые
Physical Surface(123) = {76, 70, 68, 74, 72, 110, 112, 108, 114, 116};
