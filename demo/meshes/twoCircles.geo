h=0.3;

Point(1) = {-0.1, -0.1, 0, h};
Point(2) = {-.1, 0, 0, h};
Point(3) = {0, -0.1, 0, h};
Point(4) = {-0.2, -0.1, 0, h};
Point(5) = {-0.1, -0.2, 0, h};

Circle(20) = {3, 1, 2};
Circle(21) = {2, 1, 4};
Circle(22) = {4, 1, 5};
Circle(23) = {5, 1, 3};

Point(6) = {.2, .2, 0, h};
Point(7) = {.2, .1, 0, h};
Point(8) = {.1, .2, 0, h};
Point(9) = {.2, .3, 0, h};
Point(10)= {.3, .2, 0, h};

Circle(24) = {7, 6, 8};
Circle(25) = {8, 6, 9};
Circle(26) = {9, 6, 10};
Circle(27) = {10, 6, 7};

Line Loop(100)={20,21,22,23};
Line Loop(200)={24,25,26,27};
Ruled Surface(300)={100};
Ruled Surface(400)={200};


