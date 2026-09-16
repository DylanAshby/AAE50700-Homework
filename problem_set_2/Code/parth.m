format long
O = 0;
i = 45;
w = 0;
NP = [cosd(O)*cosd(w)-sind(O)*cosd(i)*sind(w) -cosd(O)*sind(w)-sind(O)*cosd(i)*cosd(w) sind(O)*sind(i);
      sind(O)*cosd(w)+cosd(O)*cosd(i)*sind(w) -sind(O)*sind(w)+cosd(O)*cosd(i)*cosd(w) -cosd(O)*sind(i);
      sind(i)*sind(w) sind(i)*cosd(w) cosd(i)];
rhat = NP*[1;0;0]

% Velocity direction at epoch: circular orbit means the flight path angle is
% zero, so the velocity is purely transverse and lies along p-hat
vhat = NP*[0;1;0]

% Orbit-normal direction
hhat = cross(rhat,vhat)

% Apply the nadir-pointing convention: y along -rhat, x along hhat,
% and z completing the right-handed set
xhat = hhat;
yhat = -rhat;
zhat = cross(xhat,yhat);

% Stack the body axes as rows to assemble the attitude DCM at epoch
BN = [xhat'; yhat'; zhat']

% Validate the DCM using the same two tests as part (a)
orthoCheck = BN*BN'
detCheck = det(BN)

% Convert the DCM into Euler parameters (scalar last, q = [q1 q2 q3 q4])
q4 = 0.5*sqrt(trace(BN)+1);
q1 = (BN(2,3)-BN(3,2))/(4*q4);
q2 = (BN(3,1)-BN(1,3))/(4*q4);
q3 = (BN(1,2)-BN(2,1))/(4*q4);
q0 = [q1; q2; q3; q4]

% Verify the unit norm constraint and the round trip back to the DCM
%normCheck = norm(q)
%BNcheck = [1-2*q2^2-2*q3^2   2*(q1*q2+q3*q4)   2*(q1*q3-q2*q4);
%           2*(q1*q2-q3*q4)   1-2*q1^2-2*q3^2   2*(q2*q3+q1*q4);
%           2*(q1*q3+q2*q4)   2*(q2*q3-q1*q4)   1-2*q1^2-2*q2^2];
%roundTripErr = max(max(abs(BNcheck-BN)))

% Get a perturbation
wp = [3;0;0];
wp = deg2rad(wp)
tol = 10^(-12);
opt = odeset('RelTol', tol, 'AbsTol', tol);
[~, q] = ode113(@(t,q) eulerPert(q, wp), [0 1], q0, opt);

function qdot = eulerPert(q, wp)
    w1 = wp(1); w2 = wp(2); w3 = wp(3);
    qdot = 0.5*[0 w3 -w2 w1; -w3 0 w1 w2; w2 -w1 0 w3; -w1 -w2 -w3 0]*q;
end

q0=q(end,:)'
w0=[0.01;0.01;0.01];
w0=deg2rad(w0)

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

% Mean motion from the orbit parameters in Table 1
mu = 398600;                % Earth gravitational parameter [km^3/s^2]
a  = 6978;                  % semi-major axis [km]
n  = sqrt(mu/(a^3));         % [rad/s]

function xdot = rotEoM(t,x,NP,Ic,n)
    q = x(1:4);
    w = x(5:7);

    % Renormalize to hold the Euler parameter constraint against drift
    q = q/norm(q);
    w1 = w(1); w2 = w(2); w3 = w(3);

    % Current true anomaly, linear in time since the orbit is circular
    nu = n*t;

    % Rebuild the attitude DCM from the current Euler parameters
    BNq = [1-2*q(2)^2-2*q(3)^2       2*(q(1)*q(2)+q(3)*q(4))   2*(q(1)*q(3)-q(2)*q(4));
           2*(q(1)*q(2)-q(3)*q(4))   1-2*q(1)^2-2*q(3)^2       2*(q(2)*q(3)+q(1)*q(4));
           2*(q(1)*q(3)+q(2)*q(4))   2*(q(2)*q(3)-q(1)*q(4))   1-2*q(1)^2-2*q(2)^2];

    % Map the position unit vector from the P frame through N into B
    rhatB = BNq * NP * [cos(nu); sin(nu); 0];

    % Kinematics and dynamics
    qdot = 0.5*[0 w3 -w2 w1; -w3 0 w1 w2; w2 -w1 0 w3; -w1 -w2 -w3 0]*q;
    wdot = Ic\(3*(n^2)*cross(rhatB, Ic*rhatB) - cross(w, Ic*w));

    xdot = [qdot; wdot];
end

% Set LaTeX as default interpreter for all text
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultAxesTitleFontSizeMultiplier', 1.4);

% Compute the orbital period and run the simulation for several orbits
P = 2*pi/n;
tol = 10^(-12);
opt = odeset('RelTol', tol, 'AbsTol', tol, 'MaxStep', P/200);
[simtime, x] = ode113(@(t,x) rotEoM(t,x,NP,Ic,n), [0 5*P], [q0;w0], opt);

% Convert time to orbits and recover the nadir pointing error at each output
torb = simtime/P;
N = size(x,1);
err = zeros(N,1);
qnorm = zeros(N,1);
for k = 1:N
    qk = x(k,1:4)';
    qnorm(k) = norm(qk);
    qk = qk/qnorm(k);
    nu = n*simtime(k);
    BNk = [1-2*qk(2)^2-2*qk(3)^2        2*(qk(1)*qk(2)+qk(3)*qk(4))   2*(qk(1)*qk(3)-qk(2)*qk(4));
           2*(qk(1)*qk(2)-qk(3)*qk(4))  1-2*qk(1)^2-2*qk(3)^2         2*(qk(2)*qk(3)+qk(1)*qk(4));
           2*(qk(1)*qk(3)+qk(2)*qk(4))  2*(qk(2)*qk(3)-qk(1)*qk(4))   1-2*qk(1)^2-2*qk(2)^2];
    rhatBk = BNk*NP*[cos(nu); sin(nu); 0];
    err(k) = acosd(max(min(-rhatBk(2),1),-1));   % y-hat is the nadir axis
end

% Plot 1: Euler parameter time history (each component on its own subplot)
figure
subplot(4,1,1)
plot(torb, x(:,1), 'Color', [0.000 0.447 0.741])
ylabel('$q_1$ [-]')
title('Attitude (Euler parameter) time history', 'FontSize',20)
grid on
subplot(4,1,2)
plot(torb, x(:,2), 'Color', [0.850 0.325 0.098])
ylabel('$q_2$ [-]')
grid on
subplot(4,1,3)
plot(torb, x(:,3), 'Color', [0.929 0.694 0.125])
ylabel('$q_3$ [-]')
grid on
subplot(4,1,4)
plot(torb, x(:,4), 'Color', [0.494 0.184 0.556])
xlabel('Time [orbits]')
ylabel('$q_4$ [-]')
grid on

% Plot 2: body rate time history (each component on its own subplot)
figure
subplot(3,1,1)
plot(torb, x(:,5), 'Color', [0.466 0.674 0.188])
ylabel('$\omega_1$ [rad/s]')
title('Body angular rate time history',FontSize=20)
grid on
subplot(3,1,2)
plot(torb, x(:,6), 'Color', [0.301 0.745 0.933])
ylabel('$\omega_2$ [rad/s]')
grid on
subplot(3,1,3)
plot(torb, x(:,7), 'Color', [0.635 0.078 0.184])
xlabel('Time [orbits]')
ylabel('$\omega_3$ [rad/s]')
grid on

% Plot 3: deviation from nadir pointing
figure
plot(torb, err, 'Color', [0.000 0.000 0.000])
xlabel('Time [orbits]')
ylabel('Pointing error [deg]')
title('Deviation of $\hat{\mathbf{y}}$ from nadir')
grid on

% Report the bounds observed over the simulation
fprintf('Max nadir pointing error: %f deg\n', max(err))
fprintf('Max |omega| over simulation: %e rad/s\n', max(vecnorm(x(:,5:7),2,2)))
fprintf('Max |norm(q) - 1| over simulation: %e\n', max(abs(qnorm-1)))