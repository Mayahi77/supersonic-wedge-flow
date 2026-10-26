// Gmsh project created on Mon Oct 26 04:59:55 2026
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0.5, 0, 0, 1.0};
//+
Point(3) = {1.5, 0.176326980708465, 0, 1.0};
//+
Point(4) = {1.5, 1, 0, 1.0};
//+
Point(5) = {0.5, 1, 0, 1.0};
//+
Point(6) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Line(7) = {2, 5};
//+
Curve Loop(1) = {1, 7, 5, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, -4, -3, -2};
//+
Plane Surface(2) = {2};
//+
Transfinite Curve {1, 5} = 26 Using Progression 1;
//+
Transfinite Curve {2, 4} = 50 Using Progression 1;
//+
Transfinite Curve {6, 7, 3} = 50 Using Progression 1;
//+
Transfinite Surface {1} = {1, 2, 5, 6};
//+
Transfinite Surface {2} = {2, 3, 4, 5};
//+
Recombine Surface {1, 2};


//+
Extrude {0, 0, 0.01} {
  Surface{1}; Surface{2};
  Layers {1};
  Recombine;

}

//+
Physical Surface("inlet") = {6};
//+
Physical Surface("outlet") = {9};
//+
Physical Surface("lower") = {3, 10};
//+
Physical Surface("upper") = {5, 8};
//+
Physical Surface("frontAndBack") = {1, 2, 7, 11};
//+
Physical Volume("fluid") = {1, 2};

// Portable OpenFOAM import: ASCII MSH 2.2.
Mesh.MshFileVersion = 2.2;
Mesh.Binary = 0;
