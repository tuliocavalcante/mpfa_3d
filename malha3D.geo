cl1 = 1;

Point(1) = {0, 0, 0, cl1};
Point(2) = {1, 0, 0, cl1};
Point(3) = {1, 1, 0, cl1};
Point(4) = {0, 1, 0, cl1};
Point(5) = {0, 0, 1, cl1};
Point(6) = {1, 0, 1, cl1};
Point(7) = {1, 1, 1, cl1};
Point(8) = {0, 1, 1, cl1};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

Line Loop(1) = {8, 5, 6, 7};
Plane Surface(1) = {1};

Line Loop(2) = {6, -11, -2, 10};
Plane Surface(2) = {2};

Line Loop(3) = {3, 4, 1, 2};
Plane Surface(3) = {3};

Line Loop(4) = {12, 8, -9, -4};
Plane Surface(4) = {4};

Line Loop(5) = {5, -10, -1, 9};
Plane Surface(5) = {5};

Line Loop(6) = {7, -12, -3, 11};
Plane Surface(6) = {6};

Surface Loop(1) = {6, 1, 4, 5, 2, 3};
Volume(1) = {1};

Physical Point(101) = {4, 8, 5, 1};
Physical Point(102) = {3, 7, 6, 2};
Physical Line(101) = {12, 8, 4, 9};
Physical Line(102) = {11, 2, 6, 10};
Physical Line(201) = {3, 7, 1, 5};

Physical Surface(101) = {4};
Physical Surface(102) = {2};
Physical Surface(201) = {6, 3, 1, 5};
Physical Volume(1) = {1};
