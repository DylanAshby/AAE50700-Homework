%% PS9 (b.4) - Comparison of extended and non-extended phase-space solutions
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

%% Eccentric and true anomaly (M = 0 gives E = 0 and nu = 0)
E  = 0;
nu = atan2(sqrt(1 - e^2)*sin(E), cos(E) - e);

%% Orbital radius and speed
p    = a*(1 - e^2);                 % semi-latus rectum [km]
rmag = p/(1 + e*cos(nu));           % orbital radius [km]
vmag = sqrt(2*mu/rmag - mu/a);      % speed (vis-viva) [km/s]

%% Position and velocity in the RTN frame
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

%% Initial Hamiltonian and p0(0) for the extended run
R0  = dM*[cos(nM*0 + phiM); sin(nM*0 + phiM); 0];
H0  = 0.5*(p0vec.'*p0vec) - mu/norm(q0vec) - muM/norm(R0 - q0vec) + (muM/dM^3)*(q0vec.'*R0);
p00 = -H0;

%% Initial states for both formulations
x0ext = [q0vec; p0vec; 0; p00];   % extended [q; p; q0; p0]
x0reg = [q0vec; p0vec];           % non-extended [q; p]

%% Common time grid
Ngrid = 10000;
tgrid = linspace(0, TM, Ngrid);

%% Run both simulations on the common time grid
tol = 10^(-12);
opt = odeset('RelTol', tol, 'AbsTol', tol);
[~, xext] = ode45(@(t,x) threeBodyPhaseSpace(t, x, mu, muM, dM, nM, phiM), tgrid, x0ext, opt);
[~,       xreg] = ode45(@(t,x) threeBodyNonExtended(t, x, mu, muM, dM, nM, phiM), tgrid, x0reg, opt);

%% Position and velocity differences between the two solutions
dqvec = abs(xext(:,1:3) - xreg(:,1:3));   % position difference components [km]
dpvec = abs(xext(:,4:6) - xreg(:,4:6));   % velocity difference components [km/s]
dq    = vecnorm(xext(:,1:3) - xreg(:,1:3), 2, 2);   % position difference magnitude [km]
dp    = vecnorm(xext(:,4:6) - xreg(:,4:6), 2, 2);   % velocity difference magnitude [km/s]

%% Maximum differences over the common time grid
dqmax = max(dqvec);   % [max |dx|, max |dy|, max |dz|]
dpmax = max(dpvec);   % [max |dpx|, max |dpy|, max |dpz|]

%% Report the maximum differences
fprintf('Common time grid points        = %d\n\n', Ngrid);
fprintf('Max position differences:\n');
fprintf('  max |dq_x|                   = %e km\n', dqmax(1));
fprintf('  max |dq_y|                   = %e km\n', dqmax(2));
fprintf('  max |dq_z|                   = %e km\n', dqmax(3));
fprintf('  max ||dq||                   = %e km\n\n', max(dq));
fprintf('Max velocity differences:\n');
fprintf('  max |dp_x|                   = %e km/s\n', dpmax(1));
fprintf('  max |dp_y|                   = %e km/s\n', dpmax(2));
fprintf('  max |dp_z|                   = %e km/s\n', dpmax(3));
fprintf('  max ||dp||                   = %e km/s\n', max(dp));