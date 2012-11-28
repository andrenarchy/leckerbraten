// Gmsh project created on Mon Dec 26 18:27:03 2011
Point(1) = {0.1, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1.5, 0, 0.5, 1.0};
Point(4) = {1.4, 0, 0.5, 1.0};
Point(5) = {0.9, 0, 0.1, 1.0};
Point(6) = {0.1, 0, 0.1, 1.0};
Line(1) = {6, 5};
Line(2) = {5, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};
Line(5) = {2, 1};
Extrude {{0, 0, 1}, {0, 0, 0}, 0.975*Pi} {
  Line{1, 2, 3, 4, 5};
}
Extrude {{0, 0, 1}, {0, 0, 0}, -0.975*Pi} {
  Line{1, 2, 3, 4, 5};
}
Circle(46) = {40, 30, 24};
Circle(47) = {34, 12, 8};
Circle(48) = {36, 20, 14};
Circle(49) = {38, 20, 22};
Delete {
  Surface{9};
}
Delete {
  Surface{25};
}
Delete {
  Surface{29};
}
Delete {
  Surface{45};
}
Delete {
  Surface{13};
}
Delete {
  Surface{21};
}
Delete {
  Surface{17};
}
Delete {
  Surface{33};
}
Delete {
  Surface{41};
}
Delete {
  Surface{37};
}
Line Loop(52) = {40, 46, -20};
Plane Surface(53) = {52};
Line Loop(54) = {36, 49, -16};
Line Loop(55) = {12, -48, -32};
Plane Surface(56) = {54, 55};
Line Loop(57) = {10, -12, -2, 8};
Ruled Surface(58) = {57};
Line Loop(59) = {2, 32, -30, -28};
Ruled Surface(60) = {59};
Line Loop(61) = {10, -48, -30, 47};
Ruled Surface(62) = {61};
Line Loop(63) = {16, 18, -20, -4};
Ruled Surface(64) = {63};
Line Loop(65) = {18, -46, -38, 49};
Ruled Surface(66) = {65};
Line Loop(67) = {36, 38, -40, -4};
Ruled Surface(68) = {67};
Surface Loop(69) = {58, 62, 56, 68, 66, 64, 53, 60, 51};
Delete {
  Surface{66};
}
Point(43) = {-1.6, -0.1, 0.4, 1.0};
Point(44) = {-1.6, 0.1, 0.4, 1.0};
Point(45) = {-1.6, 0.1, 0.25, 1.0};
Point(46) = {-1.6, -0.1, 0.25, 1.0};
Line(70) = {38, 43};
Line(71) = {43, 46};
Line(72) = {46, 40};
Line(73) = {22, 44};
Line(74) = {44, 45};
Line(75) = {45, 24};
Line(76) = {44, 43};
Line(77) = {46, 45};
Point(47) = {-3, 0.1, 0.4, 1.0};
Point(48) = {-3, -0.1, 0.4, 1.0};
Point(49) = {-3, 0.1, 0.3, 1.0};
Point(50) = {-3, -0.1, 0.3, 1.0};
Point(51) = {-3.1, 0.0, 0.4, 1.0};
Point(52) = {-3.1, 0.0, 0.3, 1.0};
BSpline(78) = {46, 50, 52, 49, 45};
BSpline(79) = {43, 48, 51, 47, 44};
Line Loop(80) = {78, -77};
Plane Surface(81) = {80};
Line Loop(82) = {79, 76};
Plane Surface(83) = {82};
Line Loop(84) = {71, 78, -74, -79};
Ruled Surface(85) = {84};
Line Loop(86) = {38, -72, -71, -70};
Plane Surface(87) = {86};
Line Loop(88) = {18, -75, -74, -73};
Plane Surface(89) = {88};
Line Loop(90) = {73, 76, -70, 49};
Plane Surface(91) = {90};
Line Loop(92) = {75, -46, -72, 77};
Plane Surface(93) = {92};
Point(53) = {-0.2, 0.4, 0.1, 1.0};
Point(54) = {-0.1, 0.5, 0.1, 1.0};
Point(55) = {0.1, 0.6, 0.1, 1.0};
Point(56) = {0.3, 0.6, 0.1, 1.0};
Point(57) = {0.5, 0.6, 0.1, 1.0};
Point(58) = {0.6, 0.4, 0.1, 1.0};
Point(59) = {0.7, 0.3, 0.1, 1.0};
Point(60) = {0.5, 0.1, 0.1, 1.0};
Point(61) = {0.4, 0, 0.1, 1.0};
Point(62) = {0.3, -0.1, 0.1, 1.0};
Point(63) = {0.1, -0.2, 0.1, 1.0};
Point(64) = {-0.2, -0.2, 0.1, 1.0};
Point(65) = {-0.3, -0, 0.1, 1.0};
Point(66) = {-0.35, 0.3, 0.1, 1.0};
BSpline(94) = {56, 57, 58, 59, 60, 61};
BSpline(95) = {56, 55, 54, 53, 66, 65, 64, 63, 62, 61};
Line Loop(96) = {95, -94};
Plane Surface(97) = {96};
Extrude {0, 0, 0.3} {
  Surface{97};
}
Line Loop(110) = {8, -47, -28};
Plane Surface(111) = {96, 110};
Ruled Surface(112) = {65};
Surface Loop(113) = {64, 56, 68, 53, 60, 58, 62, 111, 97, 87, 93, 89, 85, 81, 83, 91};
Surface Loop(114) = {62, 58, 56, 68, 53, 64, 60, 111, 112, 97};
Volume(115) = {114};
Surface Loop(116) = {85, 87, 93, 89, 91, 83, 81, 112};
Volume(117) = {116};
Physical Volume("pan") = {115};
Physical Volume("handle") = {117};
Physical Volume("steak") = {1};
Physical Surface("surf1") = {53};
