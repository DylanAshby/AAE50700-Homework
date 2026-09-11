% define the keplerian elements at t = 0 and the grav parameter
xkep = [6978; 0.01; deg2rad(45); 0; 0; 0];
mu = 398600;

% Run the conversion
xcart = kep2cart(xkep, mu);

% Find the norm of r
r = xcart(1:3);
rnorm = norm(r);

% Define the J2 perturbation and the Earth's radius
J2 = 1.08263*(10^-3);
r0 = 6378.137;

% Compute the scalar term of J2 pert
J2scale = (-3*mu*J2*(r0^2)) / (2*(rnorm^5));

% Compute the n1 component
n1 = (1-(5*((r(3)^2)/(rnorm^2))))*r(1);

% Compute the n2 component
n2 = (1-(5*((r(3)^2)/(rnorm^2))))*r(2);

% COmpute the n3 component 
n3 = (3-(5*((r(3)^2)/(rnorm^2))))*r(3);

aj2 = J2scale*[n1; n2; n3]

% find the moon perturbation
rmoon = [384400;0;0];
mumoon = 4902.8;

amoon = -mumoon*(((r-rmoon)/(norm(r-rmoon)^3))+((rmoon)/(norm(rmoon)^3)));

% find the sun perturbation
rmoon = [1.496*(10^8);0;0];
mumoon = 1.327*(10^11);

amoon = -mumoon*(((r-rmoon)/(norm(r-rmoon)^3))+((rmoon)/(norm(rmoon)^3)));


% find the SRP at 1AU
Am = 0.01*(1/1000)^2
asrp = ((1.02*(10^14))/(149597870.66^2))*(1*(10^(-8)))*1.3


% Find drag

% Find the norm of v
v = xcart(4:6)
vnorm = norm(v)
vhat = v/vnorm
CD = 2.2;
rho = 1.454*(10^(-13))*(1000)^3


adrag = -(1/2)*CD*Am*rho*(vnorm^2)*vhat
norm(adrag)
