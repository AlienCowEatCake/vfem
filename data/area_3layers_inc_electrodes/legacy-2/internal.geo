/* Непонятная хреновина в виде паралелепипедообразного обмылка */
/*
BASE = 100000;
CENTER = {0, 0, 0};
SIZE = {6, 10, 1, 0.5};
COEFF = 0.1;
SURFACE = 100082;
PHYSICAL = 100084;
*/

CENTERX = CENTER[0];
CENTERY = CENTER[1];
CENTERZ = CENTER[2];
HEIGHT = SIZE[2];
WIDTH1 = SIZE[0];
WIDTH2 = SIZE[1];
DELTA = SIZE[3];

// Первый пошел
Point(BASE + 1) = {CENTERX + WIDTH1 / 2.0 - DELTA, CENTERY + WIDTH2 / 2.0, CENTERZ, COEFF};
Point(BASE + 2) = {CENTERX + WIDTH1 / 2.0 - DELTA, CENTERY + WIDTH2 / 2.0, CENTERZ + HEIGHT / 2.0, COEFF};
Point(BASE + 3) = {CENTERX + WIDTH1 / 2.0 - DELTA, CENTERY + WIDTH2 / 2.0 + HEIGHT / 2.0, CENTERZ, COEFF};
Circle(BASE + 1) = {BASE + 2, BASE + 1, BASE + 3};
Point(BASE + 4) = {CENTERX + WIDTH1 / 2.0 - DELTA, CENTERY + WIDTH2 / 2.0, CENTERZ - HEIGHT / 2.0, COEFF};
Circle(BASE + 2) = {BASE + 4, BASE + 1, BASE + 3};
// Второй пошел
Point(BASE + 5) = {CENTERX - WIDTH1 / 2.0 + DELTA, CENTERY + WIDTH2 / 2.0, CENTERZ, COEFF};
Point(BASE + 6) = {CENTERX - WIDTH1 / 2.0 + DELTA, CENTERY + WIDTH2 / 2.0, CENTERZ + HEIGHT / 2.0, COEFF};
Point(BASE + 7) = {CENTERX - WIDTH1 / 2.0 + DELTA, CENTERY + WIDTH2 / 2.0 + HEIGHT / 2.0, CENTERZ, COEFF};
Circle(BASE + 3) = {BASE + 6, BASE + 5, BASE + 7};
Point(BASE + 8) = {CENTERX - WIDTH1 / 2.0 + DELTA, CENTERY + WIDTH2 / 2.0, CENTERZ - HEIGHT / 2.0, COEFF};
Circle(BASE + 4) = {BASE + 8, BASE + 5, BASE + 7};
// Пришел
Line(BASE + 5) = {BASE + 2, BASE + 6};
Line(BASE + 6) = {BASE + 3, BASE + 7};
Line(BASE + 7) = {BASE + 4, BASE + 8};
Line Loop(BASE + 8) = {BASE + 5, BASE + 3, -(BASE + 6), -(BASE + 1)};
Ruled Surface(BASE + 1) = {BASE + 8};
Line Loop(BASE + 9) = {BASE + 6, -(BASE + 4), -(BASE + 7), BASE + 2};
Ruled Surface(BASE + 2) = {BASE + 9};

// Первый пошел
Point(BASE + 11) = {CENTERX + WIDTH1 / 2.0 - DELTA, CENTERY - WIDTH2 / 2.0, CENTERZ, COEFF};
Point(BASE + 12) = {CENTERX + WIDTH1 / 2.0 - DELTA, CENTERY - WIDTH2 / 2.0, CENTERZ + HEIGHT / 2.0, COEFF};
Point(BASE + 13) = {CENTERX + WIDTH1 / 2.0 - DELTA, CENTERY - WIDTH2 / 2.0 - HEIGHT / 2.0, CENTERZ, COEFF};
Circle(BASE + 11) = {BASE + 12, BASE + 11, BASE + 13};
Point(BASE + 14) = {CENTERX + WIDTH1 / 2.0 - DELTA, CENTERY - WIDTH2 / 2.0, CENTERZ - HEIGHT / 2.0, COEFF};
Circle(BASE + 12) = {BASE + 14, BASE + 11, BASE + 13};
// Второй пошел
Point(BASE + 15) = {CENTERX - WIDTH1 / 2.0 + DELTA, CENTERY - WIDTH2 / 2.0, CENTERZ, COEFF};
Point(BASE + 16) = {CENTERX - WIDTH1 / 2.0 + DELTA, CENTERY - WIDTH2 / 2.0, CENTERZ + HEIGHT / 2.0, COEFF};
Point(BASE + 17) = {CENTERX - WIDTH1 / 2.0 + DELTA, CENTERY - WIDTH2 / 2.0 - HEIGHT / 2.0, CENTERZ, COEFF};
Circle(BASE + 13) = {BASE + 16, BASE + 15, BASE + 17};
Point(BASE + 18) = {CENTERX - WIDTH1 / 2.0 + DELTA, CENTERY - WIDTH2 / 2.0, CENTERZ - HEIGHT / 2.0, COEFF};
Circle(BASE + 14) = {BASE + 18, BASE + 15, BASE + 17};
// Пришел
Line(BASE + 15) = {BASE + 12, BASE + 16};
Line(BASE + 16) = {BASE + 13, BASE + 17};
Line(BASE + 17) = {BASE + 14, BASE + 18};
Line Loop(BASE + 18) = {BASE + 15, BASE + 13, -(BASE + 16), -(BASE + 11)};
Ruled Surface(BASE + 11) = {BASE + 18};
Line Loop(BASE + 19) = {BASE + 16, -(BASE + 14), -(BASE + 17), BASE + 12};
Ruled Surface(BASE + 12) = {BASE + 19};

// Первый пошел
Point(BASE + 21) = {CENTERX + WIDTH1 / 2.0, CENTERY + WIDTH2 / 2.0 - DELTA, CENTERZ, COEFF};
Point(BASE + 22) = {CENTERX + WIDTH1 / 2.0, CENTERY + WIDTH2 / 2.0 - DELTA, CENTERZ + HEIGHT / 2.0, COEFF};
Point(BASE + 23) = {CENTERX + WIDTH1 / 2.0 + HEIGHT / 2.0, CENTERY + WIDTH2 / 2.0 - DELTA, CENTERZ, COEFF};
Circle(BASE + 21) = {BASE + 22, BASE + 21, BASE + 23};
Point(BASE + 24) = {CENTERX + WIDTH1 / 2.0, CENTERY + WIDTH2 / 2.0 - DELTA, CENTERZ - HEIGHT / 2.0, COEFF};
Circle(BASE + 22) = {BASE + 24, BASE + 21, BASE + 23};
// Второй пошел
Point(BASE + 25) = {CENTERX + WIDTH1 / 2.0, CENTERY - WIDTH2 / 2.0 + DELTA, CENTERZ, COEFF};
Point(BASE + 26) = {CENTERX + WIDTH1 / 2.0, CENTERY - WIDTH2 / 2.0 + DELTA, CENTERZ + HEIGHT / 2.0, COEFF};
Point(BASE + 27) = {CENTERX + WIDTH1 / 2.0 + HEIGHT / 2.0, CENTERY - WIDTH2 / 2.0 + DELTA, CENTERZ, COEFF};
Circle(BASE + 23) = {BASE + 26, BASE + 25, BASE + 27};
Point(BASE + 28) = {CENTERX + WIDTH1 / 2.0, CENTERY - WIDTH2 / 2.0 + DELTA, CENTERZ - HEIGHT / 2.0, COEFF};
Circle(BASE + 24) = {BASE + 28, BASE + 25, BASE + 27};
// Пришел
Line(BASE + 25) = {BASE + 22, BASE + 26};
Line(BASE + 26) = {BASE + 23, BASE + 27};
Line(BASE + 27) = {BASE + 24, BASE + 28};
Line Loop(BASE + 28) = {BASE + 25, BASE + 23, -(BASE + 26), -(BASE + 21)};
Ruled Surface(BASE + 21) = {BASE + 28};
Line Loop(BASE + 29) = {BASE + 26, -(BASE + 24), -(BASE + 27), BASE + 22};
Ruled Surface(BASE + 22) = {BASE + 29};

// Первый пошел
Point(BASE + 31) = {CENTERX - WIDTH1 / 2.0, CENTERY + WIDTH2 / 2.0 - DELTA, CENTERZ, COEFF};
Point(BASE + 32) = {CENTERX - WIDTH1 / 2.0, CENTERY + WIDTH2 / 2.0 - DELTA, CENTERZ + HEIGHT / 2.0, COEFF};
Point(BASE + 33) = {CENTERX - WIDTH1 / 2.0 - HEIGHT / 2.0, CENTERY + WIDTH2 / 2.0 - DELTA, CENTERZ, COEFF};
Circle(BASE + 31) = {BASE + 32, BASE + 31, BASE + 33};
Point(BASE + 34) = {CENTERX - WIDTH1 / 2.0, CENTERY + WIDTH2 / 2.0 - DELTA, CENTERZ - HEIGHT / 2.0, COEFF};
Circle(BASE + 32) = {BASE + 34, BASE + 31, BASE + 33};
// Второй пошел
Point(BASE + 35) = {CENTERX - WIDTH1 / 2.0, CENTERY - WIDTH2 / 2.0 + DELTA, CENTERZ, COEFF};
Point(BASE + 36) = {CENTERX - WIDTH1 / 2.0, CENTERY - WIDTH2 / 2.0 + DELTA, CENTERZ + HEIGHT / 2.0, COEFF};
Point(BASE + 37) = {CENTERX - WIDTH1 / 2.0 - HEIGHT / 2.0, CENTERY - WIDTH2 / 2.0 + DELTA, CENTERZ, COEFF};
Circle(BASE + 33) = {BASE + 36, BASE + 35, BASE + 37};
Point(BASE + 38) = {CENTERX - WIDTH1 / 2.0, CENTERY - WIDTH2 / 2.0 + DELTA, CENTERZ - HEIGHT / 2.0, COEFF};
Circle(BASE + 34) = {BASE + 38, BASE + 35, BASE + 37};
// Пришел
Line(BASE + 35) = {BASE + 32, BASE + 36};
Line(BASE + 36) = {BASE + 33, BASE + 37};
Line(BASE + 37) = {BASE + 34, BASE + 38};
Line Loop(BASE + 38) = {BASE + 35, BASE + 33, -(BASE + 36), -(BASE + 31)};
Ruled Surface(BASE + 31) = {BASE + 38};
Line Loop(BASE + 39) = {BASE + 36, -(BASE + 34), -(BASE + 37), BASE + 32};
Ruled Surface(BASE + 32) = {BASE + 39};

// Первый
Point(BASE + 41) = {CENTERX + WIDTH1 / 2.0 - DELTA, CENTERY + WIDTH2 / 2.0 - DELTA, CENTERZ, COEFF};
Circle(BASE + 41) = {BASE + 03, BASE + 41, BASE + 23};
Point(BASE + 42) = {CENTERX + WIDTH1 / 2.0 - DELTA, CENTERY + WIDTH2 / 2.0 - DELTA, CENTERZ + HEIGHT / 2.0, COEFF};
Circle(BASE + 42) = {BASE + 02, BASE + 42, BASE + 22};
Point(BASE + 43) = {CENTERX + WIDTH1 / 2.0 - DELTA, CENTERY + WIDTH2 / 2.0 - DELTA, CENTERZ - HEIGHT / 2.0, COEFF};
Circle(BASE + 43) = {BASE + 04, BASE + 43, BASE + 24};
Line Loop(BASE + 44) = {BASE + 43, BASE + 22, -(BASE + 41), -(BASE + 02)};
Ruled Surface(BASE + 45) = {BASE + 44};
Line Loop(BASE + 46) = {BASE + 42, BASE + 21, -(BASE + 41), -(BASE + 01)};
Ruled Surface(BASE + 47) = {BASE + 46};

// Второй
Point(BASE + 51) = {CENTERX + WIDTH1 / 2.0 - DELTA, CENTERY - WIDTH2 / 2.0 + DELTA, CENTERZ, COEFF};
Circle(BASE + 51) = {BASE + 13, BASE + 51, BASE + 27};
Point(BASE + 52) = {CENTERX + WIDTH1 / 2.0 - DELTA, CENTERY - WIDTH2 / 2.0 + DELTA, CENTERZ + HEIGHT / 2.0, COEFF};
Circle(BASE + 52) = {BASE + 12, BASE + 52, BASE + 26};
Point(BASE + 53) = {CENTERX + WIDTH1 / 2.0 - DELTA, CENTERY - WIDTH2 / 2.0 + DELTA, CENTERZ - HEIGHT / 2.0, COEFF};
Circle(BASE + 53) = {BASE + 14, BASE + 53, BASE + 28};
Line Loop(BASE + 54) = {BASE + 53, BASE + 24, -(BASE + 51), -(BASE + 12)};
Ruled Surface(BASE + 55) = {BASE + 54};
Line Loop(BASE + 56) = {BASE + 51, -(BASE + 23), -(BASE + 52), BASE + 11};
Ruled Surface(BASE + 57) = {BASE + 56};

// Третий
Point(BASE + 61) = {CENTERX - WIDTH1 / 2.0 + DELTA, CENTERY + WIDTH2 / 2.0 - DELTA, CENTERZ, COEFF};
Circle(BASE + 61) = {BASE + 07, BASE + 61, BASE + 33};
Point(BASE + 62) = {CENTERX - WIDTH1 / 2.0 + DELTA, CENTERY + WIDTH2 / 2.0 - DELTA, CENTERZ + HEIGHT / 2.0, COEFF};
Circle(BASE + 62) = {BASE + 06, BASE + 62, BASE + 32};
Point(BASE + 63) = {CENTERX - WIDTH1 / 2.0 + DELTA, CENTERY + WIDTH2 / 2.0 - DELTA, CENTERZ - HEIGHT / 2.0, COEFF};
Circle(BASE + 63) = {BASE + 08, BASE + 63, BASE + 34};
Line Loop(BASE + 64) = {BASE + 63, BASE + 32, -(BASE + 61), -(BASE + 04)};
Ruled Surface(BASE + 65) = {BASE + 64};
Line Loop(BASE + 66) = {BASE + 62, BASE + 31, -(BASE + 61), -(BASE + 03)};
Ruled Surface(BASE + 67) = {BASE + 66};

// Четвертый
Point(BASE + 71) = {CENTERX - WIDTH1 / 2.0 + DELTA, CENTERY - WIDTH2 / 2.0 + DELTA, CENTERZ, COEFF};
Circle(BASE + 71) = {BASE + 17, BASE + 71, BASE + 37};
Point(BASE + 72) = {CENTERX - WIDTH1 / 2.0 + DELTA, CENTERY - WIDTH2 / 2.0 + DELTA, CENTERZ + HEIGHT / 2.0, COEFF};
Circle(BASE + 72) = {BASE + 16, BASE + 72, BASE + 36};
Point(BASE + 73) = {CENTERX - WIDTH1 / 2.0 + DELTA, CENTERY - WIDTH2 / 2.0 + DELTA, CENTERZ - HEIGHT / 2.0, COEFF};
Circle(BASE + 73) = {BASE + 18, BASE + 73, BASE + 38};
Line Loop(BASE + 74) = {BASE + 73, BASE + 34, -(BASE + 71), -(BASE + 14)};
Ruled Surface(BASE + 75) = {BASE + 74};
Line Loop(BASE + 76) = {BASE + 72, BASE + 33, -(BASE + 71), -(BASE + 13)};
Ruled Surface(BASE + 77) = {BASE + 76};

// Теперь все смешиваем
Line Loop(BASE + 78) = {BASE + 17, BASE + 73, -(BASE + 37), -(BASE + 63), -(BASE + 07), BASE + 43, BASE + 27, -(BASE + 53)};
Plane Surface(BASE + 79) = {BASE + 78};
Line Loop(BASE + 80) = {BASE + 35, -(BASE + 72), -(BASE + 15), BASE + 52, -(BASE + 25), -(BASE + 42), BASE + 05, BASE + 62};
Plane Surface(BASE + 81) = {BASE + 80};
Surface Loop(SURFACE) = {BASE + 81, BASE + 31, BASE + 77, BASE + 75, BASE + 79, BASE + 12, BASE + 11, BASE + 57, BASE + 55, BASE + 22, BASE + 21, BASE + 47, BASE + 45, BASE + 02, BASE + 01, BASE + 67, BASE + 65, BASE + 32};
Volume(BASE + 83) = {SURFACE};
Physical Volume(PHYSICAL) = {BASE + 83};
