/* Некое подобие электрода */
/*
BASE = 500;
CENTER = {0, 0, 0};
SIZE = {0.5, 0.1};
COEFF = 0.1;
PHYSICAL = 1;
*/

LENGTH = SIZE[0];
RADIUS = SIZE[1];

Point(BASE + 1) = {CENTER[0], CENTER[1], CENTER[2] + (LENGTH / 2.0), COEFF};
Point(BASE + 2) = {CENTER[0] + RADIUS, CENTER[1], CENTER[2] + (LENGTH / 2.0), COEFF};
Point(BASE + 3) = {CENTER[0], CENTER[1] + RADIUS, CENTER[2] + (LENGTH / 2.0), COEFF};
Point(BASE + 4) = {CENTER[0] - RADIUS, CENTER[1], CENTER[2] + (LENGTH / 2.0), COEFF};
Point(BASE + 5) = {CENTER[0], CENTER[1] - RADIUS, CENTER[2] + (LENGTH / 2.0), COEFF};
Circle(BASE + 1) = {BASE + 2, BASE + 1, BASE + 3};
Circle(BASE + 2) = {BASE + 3, BASE + 1, BASE + 4};
Circle(BASE + 3) = {BASE + 4, BASE + 1, BASE + 5};
Circle(BASE + 4) = {BASE + 5, BASE + 1, BASE + 2};
Line(BASE + 23) = {BASE + 1, BASE + 3};
Line(BASE + 24) = {BASE + 1, BASE + 4};
Line(BASE + 25) = {BASE + 1, BASE + 5};
Line(BASE + 26) = {BASE + 1, BASE + 2};
Line Loop(BASE + 31) = {BASE + 23, BASE + 2, -(BASE + 24)};
Ruled Surface(BASE + 32) = {BASE + 31};
Line Loop(BASE + 33) = {BASE + 23, -(BASE + 1), -(BASE + 26)};
Ruled Surface(BASE + 34) = {BASE + 33};
Line Loop(BASE + 35) = {BASE + 26, -(BASE + 4), -(BASE + 25)};
Ruled Surface(BASE + 36) = {BASE + 35};
Line Loop(BASE + 37) = {BASE + 24, BASE + 3, -(BASE + 25)};
Ruled Surface(BASE + 38) = {BASE + 37};

Point(BASE + 6) = {CENTER[0], CENTER[1], CENTER[2] - (LENGTH / 2.0), COEFF};
Point(BASE + 7) = {CENTER[0] + RADIUS, CENTER[1], CENTER[2] - (LENGTH / 2.0), COEFF};
Point(BASE + 8) = {CENTER[0], CENTER[1] + RADIUS, CENTER[2] - (LENGTH / 2.0), COEFF};
Point(BASE + 9) = {CENTER[0] - RADIUS, CENTER[1], CENTER[2] - (LENGTH / 2.0), COEFF};
Point(BASE + 10) = {CENTER[0], CENTER[1] - RADIUS, CENTER[2] - (LENGTH / 2.0), COEFF};
Circle(BASE + 5) = {BASE + 7, BASE + 6, BASE + 8};
Circle(BASE + 6) = {BASE + 8, BASE + 6, BASE + 9};
Circle(BASE + 7) = {BASE + 9, BASE + 6, BASE + 10};
Circle(BASE + 8) = {BASE + 10, BASE + 6, BASE + 7};
Line(BASE + 27) = {BASE + 6, BASE + 8};
Line(BASE + 28) = {BASE + 6, BASE + 7};
Line(BASE + 29) = {BASE + 6, BASE + 10};
Line(BASE + 30) = {BASE + 6, BASE + 9};
Line Loop(BASE + 39) = {BASE + 30, BASE + 7, -(BASE + 29)};
Ruled Surface(BASE + 40) = {BASE + 39};
Line Loop(BASE + 41) = {BASE + 29, BASE + 8, -(BASE + 28)};
Ruled Surface(BASE + 42) = {BASE + 41};
Line Loop(BASE + 43) = {BASE + 30, -(BASE + 6), -(BASE + 27)};
Ruled Surface(BASE + 44) = {BASE + 43};
Line Loop(BASE + 45) = {BASE + 27, -(BASE + 5), -(BASE + 28)};
Ruled Surface(BASE + 46) = {BASE + 45};

Line(BASE + 11) = {BASE + 2, BASE + 7};
Line(BASE + 12) = {BASE + 5, BASE + 10};
Line(BASE + 13) = {BASE + 4, BASE + 9};
Line(BASE + 14) = {BASE + 3, BASE + 8};
Line Loop(BASE + 15) = {BASE + 8, -(BASE + 11), -(BASE + 4), BASE + 12};
Ruled Surface(BASE + 16) = {BASE + 15};
Line Loop(BASE + 17) = {BASE + 7, -(BASE + 12), -(BASE + 3), BASE + 13};
Ruled Surface(BASE + 18) = {BASE + 17};
Line Loop(BASE + 19) = {BASE + 6, -(BASE + 13), -(BASE + 2), BASE + 14};
Ruled Surface(BASE + 20) = {BASE + 19};
Line Loop(BASE + 21) = {BASE + 5, -(BASE + 14), -(BASE + 1), BASE + 11};
Ruled Surface(BASE + 22) = {BASE + 21};

Line(BASE + 47) = {BASE + 1, BASE + 6};
Physical Line(PHYSICAL) = {BASE + 47};

Line Loop(BASE + 48) = {BASE + 25, BASE + 12, -(BASE + 29), -(BASE + 47)};
Plane Surface(BASE + 49) = {BASE + 48};
Line Loop(BASE + 50) = {BASE + 47, BASE + 27, -(BASE + 14), -(BASE + 23)};
Plane Surface(BASE + 51) = {BASE + 50};
Line Loop(BASE + 52) = {BASE + 47, BASE + 28, -(BASE + 11), -(BASE + 26)};
Plane Surface(BASE + 53) = {BASE + 52};
Line Loop(BASE + 54) = {BASE + 30, -(BASE + 13), -(BASE + 24), BASE + 47};
Plane Surface(BASE + 55) = {BASE + 54};
Surface Loop(BASE + 56) = {BASE + 40, BASE + 18, BASE + 38, BASE + 55, BASE + 49};
Volume(BASE + 57) = {BASE + 56};
Surface Loop(BASE + 58) = {BASE + 16, BASE + 42, BASE + 36, BASE + 49, BASE + 53};
Volume(BASE + 59) = {BASE + 58};
Surface Loop(BASE + 60) = {BASE + 22, BASE + 46, BASE + 34, BASE + 51, BASE + 53};
Volume(BASE + 61) = {BASE + 60};
Surface Loop(BASE + 62) = {BASE + 20, BASE + 44, BASE + 32, BASE + 51, BASE + 55};
Volume(BASE + 63) = {BASE + 62};

