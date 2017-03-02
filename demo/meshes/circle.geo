h=0.3;

Point(1)={0.,0.,.0,h};
Point(2)={0.275,0.,0.,h};
Point(3)={0.,0.275,0.,h};
Point(4)={-0.275,0.,0.,h};
Point(5)={0.,-0.275,0.,h};

Circle(10)={2,1,3};
Circle(20)={3,1,4};
Circle(30)={4,1,5};
Circle(40)={5,1,2};

Line Loop(100)={10,20,30,40};
Ruled Surface(200)={100};

Physical Line('Dirichlet')={10,20,30,40};
Physical Surface('Circle')={200};
