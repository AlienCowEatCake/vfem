/**/

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

c1 = 12;
c2 = 15;
c3 = 30;
c4 = 80;
c41 = 60;
c5 = 120;
c6 = 700;

electrode_depth = 20;
electrode_len = 5;
electrode_cl = 1;
electrode_cube_cl = 3;

// Квадрат для электродов
Point(500) = {x0, x0, zw0, electrode_cube_cl};
Point(501) = {x0, -x0, zw0, electrode_cube_cl};
Point(502) = {-x0, x0, zw0, electrode_cube_cl};
Point(503) = {-x0, -x0, zw0, electrode_cube_cl};
Point(504) = {x0, x0, zw0 - 2*electrode_depth - electrode_len, electrode_cube_cl};
Point(505) = {x0, -x0, zw0 - 2*electrode_depth - electrode_len, electrode_cube_cl};
Point(506) = {-x0, x0, zw0 - 2*electrode_depth - electrode_len, electrode_cube_cl};
Point(507) = {-x0, -x0, zw0 - 2*electrode_depth - electrode_len, electrode_cube_cl};
Point(510) = {-x0, 0, zw0 - electrode_depth, electrode_cl};
Point(511) = {-x0, 0, zw0 - electrode_depth - electrode_len, electrode_cl};
Point(512) = {x0, 0, zw0 - electrode_depth, electrode_cl};
Point(513) = {x0, 0, zw0 - electrode_depth - electrode_len, electrode_cl};
Line(501) = {502, 500};
Line(502) = {500, 501};
Line(503) = {501, 503};
Line(504) = {503, 502};
Line(505) = {502, 506};
Line(506) = {506, 507};
Line(507) = {507, 503};
Line(508) = {507, 505};
Line(509) = {505, 501};
Line(510) = {505, 504};
Line(511) = {504, 500};
Line(512) = {504, 506};
Line(513) = {506, 511};
Line(514) = {507, 511};
Line(515) = {503, 510};
Line(516) = {502, 510};
Line(519) = {500, 512};
Line(520) = {501, 512};
Line(521) = {504, 513};
Line(522) = {505, 513};
Line(523) = {510, 511};
Line(524) = {512, 513};
Point(514) = {0, 0, zw0 - electrode_depth, electrode_cube_cl};
Point(515) = {0, 0, zw0 - electrode_depth - electrode_len, electrode_cube_cl};
Line(525) = {510, 514};
Line(526) = {514, 512};
Line(527) = {511, 515};
Line(528) = {515, 513};
Line(529) = {514, 515};
Line Loop(530) = {501, 502, 503, 504};
Line Loop(531) = {525, 529, -527, -523};
Plane Surface(532) = {531};
Line Loop(533) = {528, -524, -526, 529};
Plane Surface(534) = {533};
Line Loop(535) = {526, -520, 503, 515, 525};
Plane Surface(536) = {535};
Line Loop(537) = {528, -522, -508, 514, 527};
Plane Surface(538) = {537};
Line Loop(539) = {519, -526, -525, -516, 501};
Plane Surface(540) = {539};
Line Loop(541) = {513, 527, 528, -521, 512};
Plane Surface(542) = {541};
Line Loop(543) = {514, -523, -515, -507};
Plane Surface(544) = {543};
Line Loop(545) = {516, 523, -513, -505};
Plane Surface(546) = {545};
Line Loop(547) = {515, -516, -504};
Plane Surface(548) = {547};
Line Loop(549) = {513, -514, -506};
Plane Surface(550) = {549};
Line Loop(551) = {521, -524, -519, -511};
Plane Surface(552) = {551};
Line Loop(553) = {522, -524, -520, -509};
Plane Surface(554) = {553};
Line Loop(555) = {520, -519, 502};
Plane Surface(556) = {555};
Line Loop(557) = {521, -522, 510};
Plane Surface(558) = {557};
Plane Surface(559) = {530};
Line Loop(560) = {512, 506, 508, 510};
Plane Surface(561) = {560};
Line Loop(562) = {508, 509, 503, -507};
Plane Surface(563) = {562};
Line Loop(564) = {512, -505, 501, -511};
Plane Surface(565) = {564};
Surface Loop(566) = {548, 559, 556, 540, 536};
Volume(567) = {566};
Surface Loop(568) = {550, 561, 558, 538, 542};
Volume(569) = {568};
Surface Loop(570) = {546, 565, 552, 540, 542, 532, 534};
Volume(571) = {570};
Surface Loop(572) = {563, 554, 544, 532, 534, 536, 538};
Volume(573) = {572};

Physical Line(1) = {523};
Physical Line(2) = {524};

// Куб в воздухе вокруг петли
Point(21) = {x1, x1, za0, c2};
Point(22) = {-x1, x1, za0, c2};
Point(23) = {x1, -x1, za0, c2};
Point(24) = {-x1, -x1, za0, c2};
Line(21) = {22, 21};
Line(22) = {21, 23};
Line(23) = {23, 24};
Line(24) = {24, 22};
Line Loop(25) = {24, 21, 22, 23};
Plane Surface(26) = {530, 25};
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
Surface Loop(43) = {30, 42, 26, 40, 38, 36, 559};
Volume(44) = {43};

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
Surface Loop(89) = {84, 76, 82, 88, 86, 26, 563, 554, 544, 550, 561, 558, 546, 565, 552, 548, 556};
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
Line Loop(165) = {152, 153, 154, 151};
Plane Surface(166) = {165};
Line Loop(167) = {154, -156, -114, -157};
Plane Surface(168) = {167};
Surface Loop(169) = {160, 168, 166, 162, 164, 130, 86, 76, 82, 84, 88};
Volume(170) = {169};

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
Surface Loop(189) = {186, 184, 182, 188, 180, 166};

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
Line Loop(257) = {244, 241, 242, 243};
Plane Surface(258) = {165, 257};
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
Physical Volume(11) = {44, 142};
// Воздух - PML
Physical Volume(21) = {232};
// Вода
Physical Volume(12) = {567, 569, 571, 573, 90, 170};
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

