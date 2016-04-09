// NOSRAND_SEED=1454091950 LD_PRELOAD=./nosrand.so gmsh -3 -optimize -optimize_netgen mesh2_inc_z\=5.geo
// peter@Mercury:~/Desktop$ gmsh -version
// 2.11.0

x0 = 50;
x1 = 280;
x2 = 400;
x3 = 500;
x4 = 3000;

zg0 = 0;
zg1 = -237.5;
zg2 = -170;
zg3 = -400;
zg4 = -500;
zg5 = -3000;

za0 = 0;
za2 = 240;
za3 = 400;
za4 = 500;
za5 = 3000;

zl0 = 0;
zl1 = 2;

c00 = 3;
c01 = 2.5*Log(Fabs(zl1 - zl0)+1)+c00;
c1 = 12;
c2 = 15;
c3 = 15;
c4 = 25;
c5 = 35;
c6 = 700;

// Петля
Point(1) = {0, 0, zl1, c00};
Point(2) = {x0, 0, zl1, c00};
Point(3) = {0, x0, zl1, c00};
Point(4) = {-x0, 0, zl1, c00};
Point(5) = {0, -x0, zl1, c00};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Physical Line(1) = {1, 4, 3, 2};

// Проекция петли на границу раздела сред
Point(6) = {0, 0, zl0, c01};
Point(7) = {x0, 0, zl0, c01};
Point(8) = {0, x0, zl0, c01};
Point(9) = {-x0, 0, zl0, c01};
Point(10) = {0, -x0, zl0, c01};
Circle(5) = {7, 6, 8};
Circle(6) = {8, 6, 9};
Circle(7) = {9, 6, 10};
Circle(8) = {10, 6, 7};
Line Loop(9) = {5, 6, 7, 8};
Plane Surface(10) = {9};

// Цилиндр
Line(11) = {2, 7};
Line(12) = {5, 10};
Line(13) = {4, 9};
Line(14) = {3, 8};
Line Loop(15) = {8, -11, -4, 12};
Ruled Surface(16) = {15};
Line Loop(17) = {7, -12, -3, 13};
Ruled Surface(18) = {17};
Line Loop(19) = {6, -13, -2, 14};
Ruled Surface(20) = {19};
Line Loop(21) = {5, -14, -1, 11};
Ruled Surface(22) = {21};
Surface Loop(23) = {6, 22, 10, 20, 18, 16};
Volume(24) = {23};

// Куб в воздухе вокруг петли
Point(21) = {x1, x1, za0, c2};
Point(22) = {-x1, x1, za0, c2};
Point(23) = {x1, -x1, za0, c2};
Point(24) = {-x1, -x1, za0, c2};
Line(21) = {22, 21};
Line(22) = {21, 23};
Line(23) = {23, 24};
Line(24) = {24, 22};
Line Loop(25) = {21, 22, 23, 24};
Plane Surface(26) = {9, 25};
Point(25) = {x1, x1, za2, c3};
Point(26) = {-x1, x1, za2, c3};
Point(27) = {x1, -x1, za2, c3};
Point(28) = {-x1, -x1, za2, c3};
Line(25) = {26, 25};
Line(26) = {25, 27};
Line(27) = {27, 28};
Line(28) = {28, 26};
Line Loop(29) = {26, 27, 28, 25};
Plane Surface(30) = {29};
Line(31) = {28, 24};
Line(32) = {26, 22};
Line(33) = {21, 25};
Line(34) = {27, 23};
Line Loop(35) = {21, 33, -25, 32};
Plane Surface(36) = {35};
Line Loop(37) = {24, -32, -28, 31};
Plane Surface(38) = {37};
Line Loop(39) = {23, -31, -27, 34};
Plane Surface(40) = {39};
Line Loop(41) = {22, -34, -26, -33};
Plane Surface(42) = {41};
Surface Loop(43) = {40, 26, 38, 36, 42, 30, 6, 22, 20, 18, 16};
Volume(44) = {43};

// Куб в земле вокруг петли
Point(71) = {x1, x1, zg2, c3};
Point(72) = {-x1, x1, zg2, c3};
Point(73) = {x1, -x1, zg2, c3};
Point(74) = {-x1, -x1, zg2, c3};
Line(71) = {72, 71};
Line(72) = {71, 73};
Line(73) = {73, 74};
Line(74) = {74, 72};
Line Loop(75) = {72, 73, 74, 71};
Plane Surface(76) = {75};
Line(77) = {22, 72};
Line(78) = {74, 24};
Line(79) = {73, 23};
Line(80) = {71, 21};
Line Loop(81) = {72, 79, -22, -80};
Plane Surface(82) = {81};
Line Loop(83) = {73, 78, -23, -79};
Plane Surface(84) = {83};
Line Loop(85) = {74, -77, -24, -78};
Plane Surface(86) = {85};
Line Loop(87) = {21, -80, -71, -77};
Plane Surface(88) = {87};
Surface Loop(89) = {84, 76, 82, 88, 86, 26, 10};
Volume(90) = {89};

// Куб начала PML в воздухе
Point(111) = {x2, x2, za0, c4};
Point(112) = {-x2, x2, za0, c4};
Point(113) = {x2, -x2, za0, c4};
Point(114) = {-x2, -x2, za0, c4};
Line(111) = {112, 111};
Line(112) = {111, 113};
Line(113) = {113, 114};
Line(114) = {114, 112};
Point(121) = {x2, x2, za3, c4};
Point(122) = {-x2, x2, za3, c4};
Point(123) = {x2, -x2, za3, c4};
Point(124) = {-x2, -x2, za3, c4};
Line(121) = {122, 121};
Line(122) = {121, 123};
Line(123) = {123, 124};
Line(124) = {124, 122};
Line(125) = {123, 113};
Line(126) = {124, 114};
Line(127) = {121, 111};
Line(128) = {122, 112};
Line Loop(129) = {114, 111, 112, 113};
Plane Surface(130) = {25, 129};
Line Loop(131) = {114, -128, -124, 126};
Plane Surface(132) = {131};
Line Loop(133) = {113, -126, -123, 125};
Plane Surface(134) = {133};
Line Loop(135) = {112, -125, -122, 127};
Plane Surface(136) = {135};
Line Loop(137) = {127, -111, -128, 121};
Plane Surface(138) = {137};
Line Loop(139) = {124, 121, 122, 123};
Plane Surface(140) = {139};
Surface Loop(141) = {138, 136, 130, 134, 132, 140, 36, 42, 40, 38, 30};
Volume(142) = {141};

// Куб начала PML в земле
Point(151) = {x2, x2, zg3, c4};
Point(152) = {-x2, x2, zg3, c4};
Point(153) = {x2, -x2, zg3, c4};
Point(154) = {-x2, -x2, zg3, c4};
Line(151) = {152, 151};
Line(152) = {151, 153};
Line(153) = {153, 154};
Line(154) = {154, 152};
Line(155) = {111, 151};
Line(156) = {112, 152};
Line(157) = {154, 114};
Line(158) = {113, 153};
Line Loop(159) = {113, -157, -153, -158};
Plane Surface(160) = {159};
Line Loop(161) = {112, 158, -152, -155};
Plane Surface(162) = {161};
Line Loop(163) = {155, -151, -156, 111};
Plane Surface(164) = {163};
Line Loop(165) = {152, 153, 154, 151};
Plane Surface(166) = {165};
Line Loop(167) = {154, -156, -114, -157};
Plane Surface(168) = {167};
Surface Loop(169) = {160, 168, 166, 162, 164, 130, 86, 76, 82, 84, 88};

// Включение
BASE = 100000;
CENTER = {-100, -100, zg1};
SIZE = {400, 400, 75, 75};
COEFF = 16;
SURFACE = BASE + 82;
PHYSICAL = 14;
Merge "internal.geo";
Rotate { {1.0, 1.0, 0.0}, {CENTER[0], CENTER[1], CENTER[2]}, 5 * 3.14159 / 180.0} { Volume{BASE + 83}; }
Volume(170) = {169, SURFACE};

// Куб конца PML в воздухе
Point(201) = {x3, x3, za0, c5};
Point(202) = {-x3, x3, za0, c5};
Point(203) = {x3, -x3, za0, c5};
Point(204) = {-x3, -x3, za0, c5};
Line(201) = {202, 201};
Line(202) = {201, 203};
Line(203) = {203, 204};
Line(204) = {204, 202};
Point(211) = {x3, x3, za4, c5};
Point(212) = {-x3, x3, za4, c5};
Point(213) = {x3, -x3, za4, c5};
Point(214) = {-x3, -x3, za4, c5};
Line(211) = {212, 211};
Line(212) = {211, 213};
Line(213) = {213, 214};
Line(214) = {214, 212};
Line(215) = {201, 211};
Line(216) = {202, 212};
Line(217) = {204, 214};
Line(218) = {203, 213};
Line Loop(219) = {218, 213, -217, -203};
Plane Surface(220) = {219};
Line Loop(221) = {204, 216, -214, -217};
Plane Surface(222) = {221};
Line Loop(223) = {201, 215, -211, -216};
Plane Surface(224) = {223};
Line Loop(225) = {215, 212, -218, -202};
Plane Surface(226) = {225};
Line Loop(227) = {212, 213, 214, 211};
Plane Surface(228) = {227};
Line Loop(229) = {201, 202, 203, 204};
Plane Surface(230) = {129, 229};
Surface Loop(231) = {228, 226, 224, 230, 222, 220, 140, 132, 138, 136, 134};
Volume(232) = {231};

// Куб конца PML в земле
Point(271) = {x3, x3, zg4, c5};
Point(272) = {-x3, x3, zg4, c5};
Point(273) = {x3, -x3, zg4, c5};
Point(274) = {-x3, -x3, zg4, c5};
Line(271) = {272, 271};
Line(272) = {271, 273};
Line(273) = {273, 274};
Line(274) = {274, 272};
Line(275) = {272, 202};
Line(276) = {204, 274};
Line(277) = {273, 203};
Line(278) = {271, 201};
Line Loop(279) = {276, 274, 275, -204};
Plane Surface(280) = {279};
Line Loop(281) = {271, 278, -201, -275};
Plane Surface(282) = {281};
Line Loop(283) = {202, -277, -272, 278};
Plane Surface(284) = {283};
Line Loop(285) = {277, 203, 276, -273};
Plane Surface(286) = {285};
Line Loop(287) = {273, 274, 271, 272};
Plane Surface(288) = {287};
Surface Loop(289) = {286, 284, 288, 280, 282, 160, 168, 166, 162, 164, 230};
Volume(290) = {289};

// Воздух
Physical Volume(11) = {24, 44, 142};
// Воздух - PML
Physical Volume(21) = {232};
// Земля
Physical Volume(12) = {90, 170};
// Земля - PML
Physical Volume(22) = {290};
// Краевые
Physical Surface(51) = {286, 220, 228, 288, 284, 226, 224, 282, 280, 222};

// Бак в воздухе
Point(301) = {x4, x4, za0, c6};
Point(302) = {-x4, x4, za0, c6};
Point(303) = {x4, -x4, za0, c6};
Point(304) = {-x4, -x4, za0, c6};
Point(305) = {x4, x4, za5, c6};
Point(306) = {-x4, x4, za5, c6};
Point(307) = {x4, -x4, za5, c6};
Point(308) = {-x4, -x4, za5, c6};
Line(301) = {302, 301};
Line(302) = {301, 303};
Line(303) = {303, 304};
Line(304) = {304, 302};
Line(305) = {306, 308};
Line(306) = {308, 307};
Line(307) = {307, 305};
Line(308) = {305, 306};
Line(309) = {306, 302};
Line(310) = {301, 305};
Line(311) = {303, 307};
Line(312) = {308, 304};
Line Loop(313) = {303, -312, 306, -311};
Plane Surface(314) = {313};
Line Loop(315) = {304, -309, 305, 312};
Plane Surface(316) = {315};
Line Loop(317) = {301, 310, 308, 309};
Plane Surface(318) = {317};
Line Loop(319) = {302, 311, 307, -310};
Plane Surface(320) = {319};
Line Loop(321) = {306, 307, 308, 305};
Plane Surface(322) = {321};
Line Loop(323) = {304, 301, 302, 303};
Plane Surface(324) = {229, 323};
Surface Loop(325) = {314, 324, 320, 322, 318, 316, 228, 226, 224, 222, 220};
Volume(326) = {325};

// Бак в земле
Point(361) = {x4, x4, zg5, c6};
Point(362) = {-x4, x4, zg5, c6};
Point(363) = {x4, -x4, zg5, c6};
Point(364) = {-x4, -x4, zg5, c6};
Line(361) = {362, 361};
Line(362) = {361, 363};
Line(363) = {363, 364};
Line(364) = {364, 362};
Line(365) = {364, 304};
Line(366) = {303, 363};
Line(367) = {362, 302};
Line(368) = {301, 361};
Line Loop(369) = {365, -303, 366, 363};
Plane Surface(370) = {369};
Line Loop(371) = {365, 304, -367, -364};
Plane Surface(372) = {371};
Line Loop(373) = {367, 301, 368, -361};
Plane Surface(374) = {373};
Line Loop(375) = {362, -366, -302, 368};
Plane Surface(376) = {375};
Line Loop(377) = {364, 361, 362, 363};
Plane Surface(378) = {377};
Surface Loop(379) = {370, 372, 374, 376, 378, 324, 286, 284, 288, 280, 282};
Volume(380) = {379};

// Воздух-бак
Physical Volume(31) = {326};
// Земля-бак
Physical Volume(32) = {380};
// Краевые-бак
Physical Surface(52) = {370, 314, 316, 372, 376, 320, 318, 374, 322, 378};


