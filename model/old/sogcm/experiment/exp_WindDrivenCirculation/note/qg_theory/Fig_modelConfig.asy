import graph;

size(0,200);

real zmin=0, zmax=1;
real zoffset=0.15;
real blroffset=0.1;
real ymin=-1.3, ymax=1.3;
real yhwidth = 1.3;
pair O = (0,0);


//** Fill each region using diffrent colors

// interior
fill((ymin,zmin)--(ymin,zmax)--(ymax,zmax)--(ymax,zmin)--cycle, paleblue);
// ekman layer
fill((ymin,zmin+blroffset)--(ymin,zmax-blroffset)--(ymax,zmax-blroffset)--(ymax,zmin+blroffset)--cycle, palecyan);
// ground
fill((ymin,zmin-zoffset)--(ymax,zmin-zoffset)--(ymax,zmin)--(ymin,zmin)--cycle,lightgrey);

//** Write axes

// z-axis
yaxis("$z_*$", ymin=0, ymax=zmax+zoffset, Arrow, above=true);
draw((ymin,zmax)--(ymax,zmax));
label("$z_*=H_*$", (ymin,zmax), NE);
label("$z_*=0$", (ymin,0), SE);

// y-axis
xaxis("$y_*$", xmin=ymin, xmax=ymax, Arrow, above=true);
label("$\theta=\theta_0$", O, 3S);

// x-axis
filldraw(Circle(O,0.01*yhwidth), black, black);
draw(Circle(O,0.03*yhwidth));
label("$x_*$", O, 2NE);

//** Aditional information

// surface stress
label("$(\tau_{x*},\tau_{y*})=(\tau_{x*}(y),0)$", (yhwidth,zmax), NW);

// ekman layer
draw("Ekman layer", (-0.5*yhwidth,zmax-blroffset)--(-0.5*yhwidth,zmax),RightSide,Arrows,PenMargin);
draw("Ekman layer", (-0.5*yhwidth,zmin)--(-0.5*yhwidth,zmin+blroffset),RightSide,Arrows,PenMargin);

//interior
label("interior", (-0.5*yhwidth,0.5*zmax), E);
