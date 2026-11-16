%% PS9 (b.3) - Extended phase-space integration over one lunar revolution
clear; clc; close all;

% Set LaTeX as default interpreter for all text
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultAxesFontSize', 20);
set(groot, 'defaultTextFontSize', 20);
set(groot, 'defaultLegendFontSize', 20);

%% Constants given by the problem
mu    = 398600;      % Earth gravitational parameter [km^3/s^2]
muM   = 4902.8;      % Moon gravitational parameter [km^3/s^2]
dM    = 384400;      % Moon orbital radius [km]
phiM  = pi/4;        % Moon initial phase [rad]

%% Lunar mean motion and period
nM = sqrt((mu + muM)/dM^3);   % [rad/s]
TM = 2*pi/nM;                 % [s]

%% Given Keplerian elements of the spacecraft
a  = 42164;          % semi-major axis [km]
e  = 0.01;           % eccentricity [-]
i  = deg2rad(20);    % inclination [rad]
O  = 0;              % right ascension of the ascending node [rad]
w  = 0;              % argument of periapsis [rad]
M  = 0;              % mean anomaly [rad]

%% Eccentric and true anomaly
% M = 0 gives E = 0 directly from Kepler's equation (M = E - e*sin(E))
E  = 0;
nu = atan2(sqrt(1 - e^2)*sin(E), cos(E) - e);

%% Orbital radius and speed
p    = a*(1 - e^2);                 % semi-latus rectum [km]
rmag = p/(1 + e*cos(nu));           % orbital radius [km]
vmag = sqrt(2*mu/rmag - mu/a);      % speed (vis-viva) [km/s]

%% Position and velocity in the RTN frame (at periapsis, v is along-track)
r_rtn = [rmag; 0; 0];
v_rtn = [0; vmag; 0];

%% RTN -> inertial DCM (3-1-3 sequence with argument of latitude u)
u  = w + nu;
NB = [cos(O)*cos(u)-sin(O)*cos(i)*sin(u), -cos(O)*sin(u)-sin(O)*cos(i)*cos(u),  sin(O)*sin(i);
      sin(O)*cos(u)+cos(O)*cos(i)*sin(u), -sin(O)*sin(u)+cos(O)*cos(i)*cos(u), -cos(O)*sin(i);
      sin(i)*sin(u),                       sin(i)*cos(u),                        cos(i)];

%% Initial canonical coordinates and momenta (q = r, p = qdot = v)
q0vec = NB*r_rtn;   % q(0) [km]
p0vec = NB*v_rtn;   % p(0) [km/s]

%% Initial time coordinate
q00 = 0;

%% Moon position at q0 = 0
R0 = dM*[cos(nM*q00 + phiM); sin(nM*q00 + phiM); 0];

%% Hamiltonian at t = 0, term by term
T0      = 0.5*(p0vec.'*p0vec);                 % kinetic term
Vearth  = -mu/norm(q0vec);                     % central body term
Vmoon1  = -muM/norm(R0 - q0vec);               % direct lunar term
Vmoon2  = (muM/dM^3)*(q0vec.'*R0);             % indirect lunar term
H0      = T0 + Vearth + Vmoon1 + Vmoon2;

%% Initial momentum conjugate to time
p00 = -H0;

%% Full extended-phase-space initial state [q; p; q0; p0]
x0 = [q0vec; p0vec; q00; p00];

%% Run the simulation for one lunar revolution
tol = 10^(-12);
opt = odeset('RelTol', tol, 'AbsTol', tol, 'Refine', 10);
[simtime, x] = ode45(@(t,x) threeBodyPhaseSpace(t, x, mu, muM, dM, nM, phiM), [0 TM], x0, opt);

% Time in days for plotting
tdays = simtime/86400;

%% Compute H(t) and Hhat(t) along the trajectory
N = size(x,1);
Ht    = zeros(N,1);
Hhatt = zeros(N,1);

for k = 1:N
    qk  = x(k,1:3)';
    pk  = x(k,4:6)';
    q0k = x(k,7);
    p0k = x(k,8);

    Rk = dM*[cos(nM*q0k + phiM); sin(nM*q0k + phiM); 0];
    Ht(k)    = 0.5*(pk.'*pk) - mu/norm(qk) - muM/norm(Rk - qk) + (muM/dM^3)*(qk.'*Rk);
    Hhatt(k) = Ht(k) + p0k;
end

%% Plot 1: orbit in configuration space
figure
plot3(x(:,1), x(:,2), x(:,3))
hold on
plot3(0, 0, 0, 'k.', 'MarkerSize', 20)   % Earth center
axis equal
grid on
xlabel('$x$ [km]')
ylabel('$y$ [km]')
zlabel('$z$ [km]')
legend('Spacecraft trajectory', 'Earth center', 'Location', 'best')
title('Spacecraft trajectory over one lunar revolution')

%% Plot 2: time history of q (each component on its own subplot)
figure
subplot(3,1,1)
plot(tdays, x(:,1))
ylabel('$q_x$ [km]')
title('Position time history')
grid on
xlim([tdays(1) tdays(end)])

subplot(3,1,2)
plot(tdays, x(:,2))
ylabel('$q_y$ [km]')
grid on
xlim([tdays(1) tdays(end)])

subplot(3,1,3)
plot(tdays, x(:,3))
xlabel('Time [days]')
ylabel('$q_z$ [km]')
grid on
xlim([tdays(1) tdays(end)])

%% Plot 3: time history of p (each component on its own subplot)
figure
subplot(3,1,1)
plot(tdays, x(:,4))
ylabel('$p_x$ [km/s]')
title('Momentum time history')
grid on
xlim([tdays(1) tdays(end)])

subplot(3,1,2)
plot(tdays, x(:,5))
ylabel('$p_y$ [km/s]')
grid on
xlim([tdays(1) tdays(end)])

subplot(3,1,3)
plot(tdays, x(:,6))
xlabel('Time [days]')
ylabel('$p_z$ [km/s]')
grid on
xlim([tdays(1) tdays(end)])

%% Plot 4: H(t) - H(0)
figure
plot(tdays, Ht - Ht(1))
xlabel('Time [days]')
ylabel('$H(t)-H(0)$ [km$^2$/s$^2$]')
title('Hamiltonian variation time history')
grid on
xlim([tdays(1) tdays(end)])

%% Plot 5: Hhat(t)
figure
plot(tdays, Hhatt)
xlabel('Time [days]')
ylabel('$\widehat{H}(t)$ [km$^2$/s$^2$]')
title('Extended Hamiltonian time history')
grid on
xlim([tdays(1) tdays(end)])

%% Plot 6: Hhat(t) zoomed into the first two days
tzoom1 = 2;                      % zoom window [days]
izoom1 = tdays <= tzoom1;
figure
plot(tdays(izoom1), Hhatt(izoom1))
xlabel('Time [days]')
ylabel('$\widehat{H}(t)$ [km$^2$/s$^2$]')
title('Extended Hamiltonian time history (first two days)')
grid on
xlim([0 tzoom1])

%% Plot 7: Hhat(t) zoomed into the first six hours
tzoom2 = 0.25;                   % zoom window [days]
izoom2 = tdays <= tzoom2;
figure
plot(tdays(izoom2), Hhatt(izoom2))
xlabel('Time [days]')
ylabel('$\widehat{H}(t)$ [km$^2$/s$^2$]')
title('Extended Hamiltonian time history (first six hours)')
grid on
xlim([0 tzoom2])

%% Export figures to ../Report
outdir = fullfile('..', 'Report');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Trajectory: 3D, edge-on (X-Z and Y-Z), and top-down (X-Y) views from Plot 1
fig = figure(1); ax = fig.CurrentAxes;
view(ax, 3);
exportgraphics(fig, fullfile(outdir, 'b3_trajectory_3d.png'), 'Resolution', 300);
view(ax, 0, 0);
exportgraphics(fig, fullfile(outdir, 'b3_trajectory_edgeon_xz.png'), 'Resolution', 300);
view(ax, 90, 0);
exportgraphics(fig, fullfile(outdir, 'b3_trajectory_edgeon_yz.png'), 'Resolution', 300);
view(ax, 0, 90);
exportgraphics(fig, fullfile(outdir, 'b3_trajectory_topdown.png'), 'Resolution', 300);
view(ax, 3);

% Position and momentum time histories from Plots 2 and 3
fig = figure(2);
exportgraphics(fig, fullfile(outdir, 'b3_q.png'), 'Resolution', 300);
fig = figure(3);
exportgraphics(fig, fullfile(outdir, 'b3_p.png'), 'Resolution', 300);

% H(t) - H(0) and Hhat(t) from Plots 4 and 5
fig = figure(4);
exportgraphics(fig, fullfile(outdir, 'b3_H.png'), 'Resolution', 300);
fig = figure(5);
exportgraphics(fig, fullfile(outdir, 'b3_Hhat.png'), 'Resolution', 300);

% Zoomed Hhat(t) from Plots 6 and 7
fig = figure(6);
exportgraphics(fig, fullfile(outdir, 'b3_Hhat_zoom_2days.png'), 'Resolution', 300);
fig = figure(7);
exportgraphics(fig, fullfile(outdir, 'b3_Hhat_zoom_6hours.png'), 'Resolution', 300);

%% Print results
fprintf('n_M                 = %.10e rad/s\n', nM);
fprintf('T_M                 = %.10e s\n', TM);
fprintf('p                   = %.4f km\n', p);
fprintf('r                   = %.4f km\n', rmag);
fprintf('v                   = %.6f km/s\n', vmag);
fprintf('q(0)                = [%.4f, %.4f, %.4f] km\n', q0vec);
fprintf('p(0)                = [%.6f, %.6f, %.6f] km/s\n', p0vec);
fprintf('R(0)                = [%.4f, %.4f, %.4f] km\n', R0);
fprintf('T0                  = %.6f km^2/s^2\n', T0);
fprintf('Vearth              = %.6f km^2/s^2\n', Vearth);
fprintf('Vmoon1              = %.6f km^2/s^2\n', Vmoon1);
fprintf('Vmoon2              = %.6f km^2/s^2\n', Vmoon2);
fprintf('H(0)                = %.6f km^2/s^2\n', H0);
fprintf('p0(0)               = %.6f km^2/s^2\n', p00);
fprintf('Max |H(t) - H(0)|   = %e km^2/s^2\n', max(abs(Ht - Ht(1))));
fprintf('Max |Hhat(t)|       = %e km^2/s^2\n', max(abs(Hhatt)));
fprintf('Final q0 - T_M      = %e s\n', x(end,7) - TM);