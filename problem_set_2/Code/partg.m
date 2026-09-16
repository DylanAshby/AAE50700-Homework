% Define the parameters of the bus
m = 50;
Lx = 0.6;
Ly = 0.6;
Lz = 0.8;

% Compute the moments of inertia
Ixx = (m/12)*((Ly^2)+(Lz^2))
Iyy = (m/12)*((Lx^2)+(Lz^2))
Izz = (m/12)*((Lx^2)+(Ly^2))

% Define the bus inertia tensor
Ibus = diag([Ixx,Iyy,Izz]);

% Define the parameters of the panels
m = 5;
Lx = 0;
Ly = 1;
Lz = 1.5;

% Compute the moments of inertia
Ixx = (m/12)*((Ly^2)+(Lz^2))
Iyy = (m/12)*((Lx^2)+(Lz^2))
Izz = (m/12)*((Lx^2)+(Ly^2))

% Define the panel tensor
Ipan = diag([Ixx,Iyy,Izz]);

% Compute the offset vector
d=(0.6+1)/2*[0;1;0];

% Compute the first offset tensor
I1 = Ipan + m*((norm(d)^2)*eye(3)-(d*d'));

% Compute the second
I2 = Ipan + m*((norm(-d)^2)*eye(3)-((-d)*(-d')));

Ic=Ibus+I1+I2;