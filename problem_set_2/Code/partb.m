% Define the euler param unit vector and rotation
theta = 120;
lambda = [1;1;1]*(1/sqrt(3));

% Compute the euler parameter and verify it is unit length
euler = -[lambda*sind(theta/2); cosd(theta/2)]
norm(euler)

% Define the DCM for this euler parameter from inertial to body frame
e1 = euler(1);  e2 = euler(2);  e3 = euler(3);  e4 = euler(4);

BN = [1-2*e2^2-2*e3^2,      2*(e1*e2+e3*e4),      2*(e1*e3-e2*e4);
      2*(e1*e2-e3*e4),      1-2*e1^2-2*e3^2,      2*(e2*e3+e1*e4);
      2*(e1*e3+e2*e4),      2*(e2*e3-e1*e4),      1-2*e1^2-2*e2^2]

% Verify the matrix is a proper DCM
BN*BN'
det(BN)
