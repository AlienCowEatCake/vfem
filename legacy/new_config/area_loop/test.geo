// Gmsh project created on Sun Nov 23 09:41:41 2014

cl1 = 150;
cl2 = 50;
cl3 = 2;
Point(1) = {-1000, -1000, 0, cl2};
Point(2) = {1000, -1000, 0, cl2};
Point(3) = {1000, 1000, 0, cl2};
Point(4) = {-1000, 1000, 0, cl2};
Point(5) = {-1000, 1000, 1000, cl1};
Point(6) = {-1000, -1000, 1000, cl1};
Point(7) = {1000, -1000, 1000, cl1};
Point(8) = {1000, 1000, 1000, cl1};
Point(9) = {-1000, 1000, -1000, cl1};
Point(10) = {-1000, -1000, -1000, cl1};
Point(11) = {1000, -1000, -1000, cl1};
Point(12) = {1000, 1000, -1000, cl1};
Point(13) = {0, 0, 0, cl3};
Point(14) = {0, -50, 0, cl3};
Point(15) = {0, 50, 0, cl3};
Point(16) = {-50, 0, 0, cl3};
Point(17) = {50, 0, 0, cl3};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {4, 5};
Line(10) = {1, 6};
Line(11) = {2, 7};
Line(12) = {3, 8};
Circle(13) = {14, 13, 17};
Circle(14) = {17, 13, 15};
Circle(15) = {15, 13, 16};
Circle(16) = {16, 13, 14};
Line(17) = {1, 10};
Line(18) = {9, 10};
Line(19) = {4, 9};
Line(20) = {2, 11};
Line(21) = {11, 12};
Line(22) = {3, 12};
Line(23) = {9, 12};
Line(24) = {10, 11};
Line Loop(1) = {13, 14, 15, 16};
Plane Surface(1) = {1};
Line Loop(2) = {1, 2, 3, 4, -13, -14, -15, -16};
Plane Surface(2) = {2};
Line Loop(3) = {2, 12, -7, -11};
Plane Surface(3) = {3};
Line Loop(4) = {5, 6, 7, 8};
Plane Surface(4) = {4};
Line Loop(5) = {4, 10, -5, -9};
Plane Surface(5) = {5};
Line Loop(6) = {1, 11, -6, -10};
Plane Surface(6) = {6};
Line Loop(7) = {3, 9, -8, -12};
Plane Surface(7) = {7};

// Surface Loop(1) = {1, 2, 3, 4, 5, 6, 7};

Line Loop(8) = {4, 17, -18, -19};
Plane Surface(8) = {8};
Line Loop(9) = {21, -23, 18, 24};
Plane Surface(9) = {9};
Line Loop(10) = {21, -22, -2, 20};
Plane Surface(10) = {10};
Line Loop(11) = {20, -24, -17, 1};
Plane Surface(11) = {11};
Line Loop(12) = {3, 19, 23, -22};
Plane Surface(12) = {12};

Surface Loop(1) = {1, 2, 3, 4, 5, 6, 7};
Surface Loop(2) = {1, 2, 8, 9, 10, 11, 12};
Volume(1) = {1};
Volume(2) = {2};

Physical Line(1) = {13, 14, 15, 16};
Physical Surface(1) = {3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
Physical Volume(1) = {1, 2};






