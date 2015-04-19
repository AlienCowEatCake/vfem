/*
BASE = 2200;
PREV = 1200;
MULT = 3;
COEF = k * k * k;
PHYS = 40;
*/

// Куб конца PML в воздухе
Point(BASE + 1) = {(x2 + MULT * d), (x2 + MULT * d), za0, c4 * COEF};
Point(BASE + 2) = {-(x2 + MULT * d), (x2 + MULT * d), za0, c4 * COEF};
Point(BASE + 3) = {(x2 + MULT * d), -(x2 + MULT * d), za0, c4 * COEF};
Point(BASE + 4) = {-(x2 + MULT * d), -(x2 + MULT * d), za0, c4 * COEF};
Line(BASE + 1) = {BASE + 2, BASE + 1};
Line(BASE + 2) = {BASE + 1, BASE + 3};
Line(BASE + 3) = {BASE + 3, BASE + 4};
Line(BASE + 4) = {BASE + 4, BASE + 2};
Point(BASE + 11) = {(x2 + MULT * d), (x2 + MULT * d), za_p + MULT * d, c4 * COEF};
Point(BASE + 12) = {-(x2 + MULT * d), (x2 + MULT * d), za_p + MULT * d, c4 * COEF};
Point(BASE + 13) = {(x2 + MULT * d), -(x2 + MULT * d), za_p + MULT * d, c4 * COEF};
Point(BASE + 14) = {-(x2 + MULT * d), -(x2 + MULT * d), za_p + MULT * d, c4 * COEF};
Line(BASE + 11) = {BASE + 12, BASE + 11};
Line(BASE + 12) = {BASE + 11, BASE + 13};
Line(BASE + 13) = {BASE + 13, BASE + 14};
Line(BASE + 14) = {BASE + 14, BASE + 12};
Line(BASE + 15) = {BASE + 1, BASE + 11};
Line(BASE + 16) = {BASE + 2, BASE + 12};
Line(BASE + 17) = {BASE + 4, BASE + 14};
Line(BASE + 18) = {BASE + 3, BASE + 13};
Line Loop(BASE + 19) = {BASE + 18, BASE + 13, -(BASE + 17), -(BASE + 3)};
Plane Surface(BASE + 20) = {BASE + 19};
Line Loop(BASE + 21) = {BASE + 4, BASE + 16, -(BASE + 14), -(BASE + 17)};
Plane Surface(BASE + 22) = {BASE + 21};
Line Loop(BASE + 23) = {BASE + 1, BASE + 15, -(BASE + 11), -(BASE + 16)};
Plane Surface(BASE + 24) = {BASE + 23};
Line Loop(BASE + 25) = {BASE + 15, BASE + 12, -(BASE + 18), -(BASE + 2)};
Plane Surface(BASE + 26) = {BASE + 25};
Line Loop(BASE + 27) = {BASE + 12, BASE + 13, BASE + 14, BASE + 11};
Plane Surface(BASE + 28) = {BASE + 27};
Line Loop(BASE + 29) = {BASE + 1, BASE + 2, BASE + 3, BASE + 4};
Plane Surface(BASE + 30) = {PREV + 29, BASE + 29};
Surface Loop(BASE + 31) = {BASE + 28, BASE + 26, BASE + 24, BASE + 30, BASE + 22, BASE + 20, PREV + 24, PREV + 26, PREV + 28, PREV + 20, PREV + 22};
Volume(BASE + 32) = {BASE + 31};

// Куб конца PML в воде
Point(BASE + 41) = {(x2 + MULT * d), (x2 + MULT * d), zw_p - MULT * d, c4 * COEF};
Point(BASE + 42) = {-(x2 + MULT * d), (x2 + MULT * d), zw_p - MULT * d, c4 * COEF};
Point(BASE + 43) = {(x2 + MULT * d), -(x2 + MULT * d), zw_p - MULT * d, c4 * COEF};
Point(BASE + 44) = {-(x2 + MULT * d), -(x2 + MULT * d), zw_p - MULT * d, c4 * COEF};
Line(BASE + 41) = {BASE + 42, BASE + 41};
Line(BASE + 42) = {BASE + 41, BASE + 43};
Line(BASE + 43) = {BASE + 43, BASE + 44};
Line(BASE + 44) = {BASE + 44, BASE + 42};
Line(BASE + 45) = {BASE + 41, BASE + 1};
Line(BASE + 46) = {BASE + 42, BASE + 2};
Line(BASE + 47) = {BASE + 44, BASE + 4};
Line(BASE + 48) = {BASE + 43, BASE + 3};
Line Loop(BASE + 49) = {BASE + 46, BASE + 1, -(BASE + 45), -(BASE + 41)};
Plane Surface(BASE + 50) = {BASE + 49};
Line Loop(BASE + 51) = {BASE + 45, BASE + 2, -(BASE + 48), -(BASE + 42)};
Plane Surface(BASE + 52) = {BASE + 51};
Line Loop(BASE + 53) = {BASE + 48, BASE + 3, -(BASE + 47), -(BASE + 43)};
Plane Surface(BASE + 54) = {BASE + 53};
Line Loop(BASE + 55) = {BASE + 47, BASE + 4, -(BASE + 46), -(BASE + 44)};
Plane Surface(BASE + 56) = {BASE + 55};
Line Loop(BASE + 57) = {BASE + 44, BASE + 41, BASE + 42, BASE + 43};
Plane Surface(BASE + 58) = {BASE + 57};
Surface Loop(BASE + 59) = {BASE + 58, BASE + 56, BASE + 54, BASE + 52, BASE + 50, PREV + 58, PREV + 56, PREV + 54, PREV + 52, PREV + 50, BASE + 30};
Volume(BASE + 60) = {BASE + 59};

// Воздух - PML
Physical Volume(PHYS + 1) = {BASE + 32};
// Вода - PML
Physical Volume(PHYS + 2) = {BASE + 60};
// Краевые
Physical Surface(PHYS + 5) = {BASE + 24, BASE + 50, BASE + 22, BASE + 56, BASE + 54, BASE + 20, BASE + 26, BASE + 52, BASE + 58, BASE + 28};

