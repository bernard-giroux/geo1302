
lc = 0.4;
rad = 1.0;

// Sphere
Point(57) = {0.0,0.0,0.0,lc};
Point(58) = {rad,0.0,0.0,lc};
Point(59) = {0,rad,0.0,lc};
Point(60) = {-rad,0,0.0,lc};
Point(61) = {0,-rad,0.0,lc};
Point(62) = {0,0,-rad,lc};
Point(63) = {0,0,rad,lc};

Circle(108) = {58,57,59};
Circle(109) = {59,57,60};
Circle(110) = {60,57,61};
Circle(111) = {61,57,58};
Circle(112) = {59,57,62};
Circle(113) = {62,57,61};
Circle(114) = {61,57,63};
Circle(115) = {63,57,59};
Circle(116) = {58,57,63};
Circle(117) = {63,57,60};
Circle(118) = {60,57,62};
Circle(119) = {62,57,58};


Line Loop(402) = { 109, 115,-117};
Line Loop(403) = { 117, 110, 114};
Line Loop(404) = {-115,-116, 108};
Line Loop(405) = {-118,-109, 112};
Line Loop(406) = {-112,-119,-108};
Line Loop(407) = {-110, 118, 113};
Line Loop(408) = {-114, 111, 116};
Line Loop(409) = {-111, 119,-113};


Ruled Surface(602) = {402};
Ruled Surface(603) = {403};
Ruled Surface(604) = {404};
Ruled Surface(605) = {405};
Ruled Surface(606) = {406};
Ruled Surface(607) = {407};
Ruled Surface(608) = {408};
Ruled Surface(609) = {409};



Surface Loop(700) = {602,603,604,605,606,607,608,609};


Volume(900) = {700};  // Sphere

Physical Volume("Sphere") = {900};
