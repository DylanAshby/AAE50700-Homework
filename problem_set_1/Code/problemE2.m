% Define the mean anomaly and eccentricity
M = deg2rad(60);
e = 0.4;

% Define a function handle for the kepler equation and its gradient
f = @(E) E - e*sin(E) - M;
namblaf = @(E) 1 - e*cos(E);

% Define vectors for the iteration values to be stored in
E = zeros(1,10001);
fE = zeros(1,10001);
namblafE = zeros(1,10001);
deltaE = zeros(1,10001);

% Define some initial values
E(1) = M;              % Set the initial guess as mean anomaly
targetDE = 10^(-8);    % Set the target change in iteration

% Loop until change in iteration meets target
for n = 1:10000

    % Calculate this iteration of Kepler's equation and its gradient
    fE(n) = f(E(n));
    namblafE(n) = namblaf(E(n));

    % Find the next iteration of eccentric anomaly
    E(n+1) = E(n) - fE(n)/namblafE(n);

    % Calculate the difference in iteration
    deltaE(n) = abs(E(n+1) - E(n));

    % Break the loop if the difference meets the target
    if (deltaE(n) <= targetDE)
        break
    end
end