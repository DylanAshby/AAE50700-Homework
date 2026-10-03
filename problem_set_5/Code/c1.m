% Define the relevant values
a = 6628;
e = 0.01;
mu = 398600;

% Compute r and v magnitude at periapsis
p = a*(1-e^2)
r = p/(1+(e*cos(0)))
v = sqrt(((2*mu)/r)-(mu/a))

% Define the vector forms in RTN
rB = [r; 0; 0];
vB = [0; v; 0];

% Define the DCM
NB = [1 0 0; 0 cosd(45) -sind(45); 0 sind(45) cosd(45)]

% Then move rB and vB to ECI frame
rN = NB*rB
vN = NB*vB

% Compute spherical initial conditions
x = rN(1); y = rN(2); z = rN(3); xd = vN(1); yd = vN(2); zd = vN(3);
rS = [r; atan2d(y,x); atan2d(z,sqrt(x^2+y^2))]
vS = [((x*xd)+(y*yd)+(z*zd))/r;
      rad2deg(((x*yd)-(y*xd))/((x^2)+(y^2)));
      rad2deg(((((x^2)+(y^2))*zd)-(z*((x*xd)+(y*yd))))/(sqrt((x^2)+(y^2))*(r^2)))]

