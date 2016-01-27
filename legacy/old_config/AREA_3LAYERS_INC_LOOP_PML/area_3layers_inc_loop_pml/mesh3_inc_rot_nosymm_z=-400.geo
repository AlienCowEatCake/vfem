x0 = 50;
x1 = 280;
x2 = 600;
x3 = 700;
x4 = 3000;

zw0 = 0;
zw1 = -300;
zw2 = -600;

zg0 = -600;
zg1 = -720;
zg2 = -900;
zg3 = -1000;
zg4 = -3000;

za0 = 0;
za2 = 280;
za3 = 600;
za4 = 700;
za5 = 3000;

zl0 = -600;
zl1 = -400;

c00 = 3;
c01 = 2.48534*Log(Fabs(zl1 - zl0));
c1 = 12;
c2 = 15;
c3 = 30;
c4 = 80;
c41 = 60;
c5 = 120;
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

// Куб в воде вокруг петли
Point(21) = {x1, x1, zw1, c2};
Point(22) = {-x1, x1, zw1, c2};
Point(23) = {x1, -x1, zw1, c2};
Point(24) = {-x1, -x1, zw1, c2};
Line(21) = {22, 21};
Line(22) = {21, 23};
Line(23) = {23, 24};
Line(24) = {24, 22};
Point(25) = {x1, x1, zw2, c3};
Point(26) = {-x1, x1, zw2, c3};
Point(27) = {x1, -x1, zw2, c3};
Point(28) = {-x1, -x1, zw2, c3};
Line(25) = {26, 25};
Line(26) = {25, 27};
Line(27) = {27, 28};
Line(28) = {28, 26};
Line(29) = {28, 24};
Line(30) = {22, 26};
Line(31) = {25, 21};
Line(32) = {27, 23};
Line Loop(33) = {21, 22, 23, 24};
Plane Surface(34) = {33};
Line Loop(35) = {31, -21, 30, 25};
Plane Surface(36) = {35};
Line Loop(37) = {22, -32, -26, 31};
Plane Surface(38) = {37};
Line Loop(39) = {32, 23, -29, -27};
Plane Surface(40) = {39};
Line Loop(41) = {29, 24, 30, -28};
Plane Surface(42) = {41};
Line Loop(43) = {25, 26, 27, 28};
Plane Surface(44) = {9, 43};
Surface Loop(45) = {36, 38, 34, 40, 42, 44, 6, 22, 20, 18, 16};
Volume(46) = {45};

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
Plane Surface(130) = {129};
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
Surface Loop(141) = {138, 136, 130, 132, 140, 134};
Volume(142) = {141};

// Куб начала PML в воде
Point(151) = {x2, x2, zw2, c41};
Point(152) = {-x2, x2, zw2, c41};
Point(153) = {x2, -x2, zw2, c41};
Point(154) = {-x2, -x2, zw2, c41};
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
Line Loop(167) = {154, -156, -114, -157};
Plane Surface(168) = {167};
Line Loop(169) = {152, 153, 154, 151};
Plane Surface(170) = {43, 169};
Surface Loop(171) = {164, 162, 160, 168, 170, 38, 34, 36, 42, 40, 130};
Volume(172) = {171};

// Куб начала PML в грунте
Point(171) = {x2, x2, zg2, c4};
Point(172) = {-x2, x2, zg2, c4};
Point(173) = {x2, -x2, zg2, c4};
Point(174) = {-x2, -x2, zg2, c4};
Line(171) = {172, 171};
Line(172) = {171, 173};
Line(173) = {173, 174};
Line(174) = {174, 172};
Line(175) = {151, 171};
Line(176) = {152, 172};
Line(177) = {174, 154};
Line(178) = {153, 173};
Line Loop(179) = {176, -174, 177, 154};
Plane Surface(180) = {179};
Line Loop(181) = {171, -175, -151, 176};
Plane Surface(182) = {181};
Line Loop(183) = {175, 172, -178, -152};
Plane Surface(184) = {183};
Line Loop(185) = {178, 173, 177, -153};
Plane Surface(186) = {185};
Line Loop(187) = {171, 172, 173, 174};
Plane Surface(188) = {187};
Surface Loop(189) = {188, 182, 184, 186, 180, 170, 44, 10};

// Нефть!
BASE = 100000;
CENTER = {0, 0, zg1};
SIZE = {400, 400, 75, 75};
COEFF = 25;
SURFACE = BASE + 82;
PHYSICAL = 14;
Merge "internal.geo";
Rotate { {1.0, 1.0, 0.0}, {CENTER[0], CENTER[1], CENTER[2]}, 5 * 3.14159 / 180.0} { Volume{BASE + 83}; }
Volume(190) = {189, SURFACE};

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

// Куб конца PML в воде
Point(241) = {x3, x3, zw2, c5};
Point(242) = {-x3, x3, zw2, c5};
Point(243) = {x3, -x3, zw2, c5};
Point(244) = {-x3, -x3, zw2, c5};
Line(241) = {242, 241};
Line(242) = {241, 243};
Line(243) = {243, 244};
Line(244) = {244, 242};
Line(245) = {241, 201};
Line(246) = {242, 202};
Line(247) = {244, 204};
Line(248) = {243, 203};
Line Loop(249) = {246, 201, -245, -241};
Plane Surface(250) = {249};
Line Loop(251) = {245, 202, -248, -242};
Plane Surface(252) = {251};
Line Loop(253) = {248, 203, -247, -243};
Plane Surface(254) = {253};
Line Loop(255) = {247, 204, -246, -244};
Plane Surface(256) = {255};
Line Loop(257) = {242, 243, 244, 241};
Plane Surface(258) = {169, 257};
Surface Loop(259) = {256, 254, 252, 250, 258, 230, 168, 164, 162, 160};
Volume(260) = {259};

// Куб конца PML в грунте
Point(271) = {x3, x3, zg3, c5};
Point(272) = {-x3, x3, zg3, c5};
Point(273) = {x3, -x3, zg3, c5};
Point(274) = {-x3, -x3, zg3, c5};
Line(271) = {272, 271};
Line(272) = {271, 273};
Line(273) = {273, 274};
Line(274) = {274, 272};
Line(275) = {271, 241};
Line(276) = {272, 242};
Line(277) = {274, 244};
Line(278) = {273, 243};
Line Loop(279) = {276, 241, -275, -271};
Plane Surface(280) = {279};
Line Loop(281) = {275, 242, -278, -272};
Plane Surface(282) = {281};
Line Loop(283) = {278, 243, -277, -273};
Plane Surface(284) = {283};
Line Loop(285) = {277, 244, -276, -274};
Plane Surface(286) = {285};
Line Loop(287) = {274, 271, 272, 273};
Plane Surface(288) = {287};
Surface Loop(289) = {286, 284, 282, 280, 288, 180, 182, 188, 184, 186, 258};
Volume(290) = {289};

// Воздух
Physical Volume(11) = {142};
// Воздух - PML
Physical Volume(21) = {232};
// Вода
Physical Volume(12) = {24, 46, 172};
// Вода - PML
Physical Volume(22) = {260};
// Грунт
Physical Volume(13) = {190};
// Грунт - PML
Physical Volume(23) = {290};
// У нефти задано 14, см выше
// Краевые
Physical Surface(51) = {228, 222, 256, 286, 288, 282, 252, 226, 224, 220, 250, 280, 284, 254};

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

// Бак в воде
Point(331) = {x4, x4, zw2, c6};
Point(332) = {-x4, x4, zw2, c6};
Point(333) = {x4, -x4, zw2, c6};
Point(334) = {-x4, -x4, zw2, c6};
Line(331) = {332, 331};
Line(332) = {331, 333};
Line(333) = {333, 334};
Line(334) = {334, 332};
Line(335) = {333, 303};
Line(336) = {334, 304};
Line(337) = {301, 331};
Line(338) = {332, 302};
Line Loop(339) = {301, 337, -331, 338};
Plane Surface(340) = {339};
Line Loop(341) = {338, -304, -336, 334};
Plane Surface(342) = {341};
Line Loop(343) = {333, 336, -303, -335};
Plane Surface(344) = {343};
Line Loop(345) = {302, -335, -332, -337};
Plane Surface(346) = {345};
Line Loop(347) = {334, 331, 332, 333};
Plane Surface(348) = {257, 347};
Surface Loop(349) = {342, 340, 346, 344, 348, 324, 256, 254, 252, 250};
Volume(350) = {349};

// Бак в грунте
Point(361) = {x4, x4, zg4, c6};
Point(362) = {-x4, x4, zg4, c6};
Point(363) = {x4, -x4, zg4, c6};
Point(364) = {-x4, -x4, zg4, c6};
Line(361) = {362, 361};
Line(362) = {361, 363};
Line(363) = {363, 364};
Line(364) = {364, 362};
Line(365) = {363, 333};
Line(366) = {364, 334};
Line(367) = {331, 361};
Line(368) = {362, 332};
Line Loop(369) = {331, 367, -361, 368};
Plane Surface(370) = {369};
Line Loop(371) = {368, -334, -366, 364};
Plane Surface(372) = {371};
Line Loop(373) = {363, 366, -333, -365};
Plane Surface(374) = {373};
Line Loop(375) = {332, -365, -362, -367};
Plane Surface(376) = {375};
Line Loop(377) = {362, 363, 364, 361};
Plane Surface(378) = {377};
Surface Loop(379) = {372, 370, 376, 374, 378, 348, 286, 284, 282, 280, 288};
Volume(380) = {379};

// Воздух-бак
Physical Volume(31) = {326};
// Вода-бак
Physical Volume(32) = {350};
// Грунт-бак
Physical Volume(33) = {380};
// Краевые-бак
Physical Surface(52) = {322, 378, 372, 342, 316, 320, 344, 346, 376, 374, 314, 318, 340, 370};
