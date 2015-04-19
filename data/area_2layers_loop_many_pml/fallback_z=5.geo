x0 = 50;
x1 = 280;
x2 = 500;
x4 = 1500;

zw0 = 0;
zw1 = -280;
zw_p = -500;

za0 = 0;
za1 = 5;
za2 = 280;
za_p = 500;

c0 = 5;
c1 = 25;
c2 = 30;
c3 = 50;
c4 = 150;

/*
k = 1.07;
d = 50;
*/
/**/
k = 1.13;
d = 100;
/**/

// Петля
Point(1) = {0, 0, za1, c0};
Point(2) = {x0, 0, za1, c0};
Point(3) = {0, x0, za1, c0};
Point(4) = {-x0, 0, za1, c0};
Point(5) = {0, -x0, za1, c0};
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Physical Line(1) = {1, 4, 3, 2};

// Квадрат вокруг петли
Point(11) = {x1, x1, za1, c1};
Point(12) = {-x1, x1, za1, c1};
Point(13) = {x1, -x1, za1, c1};
Point(14) = {-x1, -x1, za1, c1};
Line(11) = {12, 11};
Line(12) = {11, 13};
Line(13) = {13, 14};
Line(14) = {14, 12};
Line Loop(15) = {12, 13, 14, 11};
Plane Surface(16) = {5, 15};

// Куб в воздухе вокруг петли
Point(21) = {x1, x1, za0, c2};
Point(22) = {-x1, x1, za0, c2};
Point(23) = {x1, -x1, za0, c2};
Point(24) = {-x1, -x1, za0, c2};
Line(21) = {22, 21};
Line(22) = {21, 23};
Line(23) = {23, 24};
Line(24) = {24, 22};
Line Loop(25) = {22, 23, 24, 21};
Plane Surface(26) = {25};
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
Line(31) = {21, 11};
Line(32) = {11, 25};
Line(33) = {22, 12};
Line(34) = {12, 26};
Line(35) = {23, 13};
Line(36) = {13, 27};
Line(37) = {28, 14};
Line(38) = {14, 24};
Line Loop(39) = {21, 31, -11, -33};
Plane Surface(40) = {39};
Line Loop(41) = {11, 32, -25, -34};
Plane Surface(42) = {41};
Line Loop(43) = {24, 33, -14, 38};
Plane Surface(44) = {43};
Line Loop(45) = {37, 14, 34, -28};
Plane Surface(46) = {45};
Line Loop(47) = {38, -23, 35, 13};
Plane Surface(48) = {47};
Line Loop(49) = {36, 27, 37, -13};
Plane Surface(50) = {49};
Line Loop(51) = {35, -12, -31, 22};
Plane Surface(52) = {51};
Line Loop(53) = {32, 26, -36, -12};
Plane Surface(54) = {53};
Surface Loop(55) = {30, 54, 42, 46, 50, 16, 6};
Volume(56) = {55};
Surface Loop(57) = {48, 44, 26, 52, 40, 16, 6};
Volume(58) = {57};

// Куб в воде вокруг петли
Point(71) = {x1, x1, zw1, c3};
Point(72) = {-x1, x1, zw1, c3};
Point(73) = {x1, -x1, zw1, c3};
Point(74) = {-x1, -x1, zw1, c3};
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
Surface Loop(89) = {86, 76, 82, 84, 88, 26};
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
Point(121) = {x2, x2, za_p, c4};
Point(122) = {-x2, x2, za_p, c4};
Point(123) = {x2, -x2, za_p, c4};
Point(124) = {-x2, -x2, za_p, c4};
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
Surface Loop(141) = {140, 132, 130, 134, 136, 138, 30, 54, 42, 46, 50, 40, 52, 48, 44};
Volume(142) = {141};

// Куб начала PML в воде
Point(151) = {x2, x2, zw_p, c4};
Point(152) = {-x2, x2, zw_p, c4};
Point(153) = {x2, -x2, zw_p, c4};
Point(154) = {-x2, -x2, zw_p, c4};
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
Surface Loop(169) = {166, 162, 160, 168, 164, 130, 88, 82, 76, 84, 86};
Volume(170) = {169};

// Воздух
Physical Volume(11) = {58, 56, 142};
// Вода
Physical Volume(12) = {90, 170};

// ====================================================================================================================

// Куб конца PML в воздухе
Point(201) = {(x2 + d), (x2 + d), za0, c4 * k};
Point(202) = {-(x2 + d), (x2 + d), za0, c4 * k};
Point(203) = {(x2 + d), -(x2 + d), za0, c4 * k};
Point(204) = {-(x2 + d), -(x2 + d), za0, c4 * k};
Line(201) = {202, 201};
Line(202) = {201, 203};
Line(203) = {203, 204};
Line(204) = {204, 202};
Point(211) = {(x2 + d), (x2 + d), za_p + d, c4 * k};
Point(212) = {-(x2 + d), (x2 + d), za_p + d, c4 * k};
Point(213) = {(x2 + d), -(x2 + d), za_p + d, c4 * k};
Point(214) = {-(x2 + d), -(x2 + d), za_p + d, c4 * k};
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
Point(241) = {(x2 + d), (x2 + d), zw_p - d, c4 * k};
Point(242) = {-(x2 + d), (x2 + d), zw_p - d, c4 * k};
Point(243) = {(x2 + d), -(x2 + d), zw_p - d, c4 * k};
Point(244) = {-(x2 + d), -(x2 + d), zw_p - d, c4 * k};
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
Line Loop(257) = {244, 241, 242, 243};
Plane Surface(258) = {257};
Surface Loop(259) = {258, 256, 254, 252, 250, 166, 162, 160, 168, 164, 230};
Volume(260) = {259};

// Воздух - PML
Physical Volume(21) = {232};
// Вода - PML
Physical Volume(22) = {260};
// Краевые
Physical Surface(25) = {224, 250, 222, 256, 254, 220, 226, 252, 258, 228};

// ====================================================================================================================

/*
BASE = 1200;
PREV = 200;
MULT = 2;
COEF = k * k;
PHYS = 30;
Merge "layer.geo";
*/

NUMB = Round((x4 - x2) / d) - 1;
For i In {1:NUMB:1}
  BASE = i * 1000 + 200;
  PREV = (i - 1) * 1000 + 200;
  MULT = i + 1;
  COEF = k ^ MULT;
  PHYS = (i + 2) * 10;
  Merge "layer.geo";
EndFor

Physical Surface(664828) = {BASE + 24, BASE + 50, BASE + 22, BASE + 56, BASE + 54, BASE + 20, BASE + 26, BASE + 52, BASE + 58, BASE + 28};

