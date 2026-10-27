%% Simulation

% Define the beta angle functions
bm   = deg2rad(60);
tm   = 60;
bt   = @(t)  bm*(10*(t/tm).^3-15*(t/tm).^4+6*(t/tm).^5).*(t<=tm)+bm*(t>tm);
bdt  = @(t) (bm/tm)*(30*(t/tm).^2-60*(t/tm).^3 + 30*(t/tm).^4).*(t<=tm);
bddt = @(t) (bm/tm^2)*(60*(t/tm)-180*(t/tm).^2+120*(t/tm).^3).*(t<=tm);

% Define the system constant parameters
IB  = diag([80  100  60]);
IA0 = diag([6   30   35]);
Jw  = diag([0.1 0.12 0.08]);
mA  = 24;
d   = 1.5;

% Define the initial condition of the simulation and a constant zero torque
x0 = [pi/3; pi/3; pi/3; zeros(6,1); 300; -200; 100];
u = zeros(3,1);

% Run the simulation
tspan = [0 100];
tol = 10^(-12);
opt = odeset('RelTol', tol, 'AbsTol', tol);
[simtime, xsim] = ode45(@(t,x) rwEoM(x, t, u, bt, bdt, bddt, IB, IA0, ...
                                     Jw, mA, d), tspan, x0, opt);


%% Plotting

% Set LaTeX as default interpreter for all text
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultAxesFontSize', 20);
set(groot, 'defaultTextFontSize', 20);
set(groot, 'defaultLegendFontSize', 20);

% Euler angle time history
figure
subplot(3,1,1)
plot(simtime, rad2deg(xsim(:,1)))
ylabel('$\theta_1$ [deg]')
title('Euler Angle Time History Zero Motor Torque')
grid on
subplot(3,1,2)
plot(simtime, rad2deg(xsim(:,2)))
ylabel('$\theta_2$ [deg]')
grid on
subplot(3,1,3)
plot(simtime, rad2deg(xsim(:,3)))
xlabel('Time [s]')
ylabel('$\theta_3$ [deg]')
grid on

% RW angular rate time history
figure
subplot(3,1,1)
plot(simtime, xsim(:,10))
ylabel('$\Omega_1$ [rad/s]')
title('RW Angular Rate Time History Zero Motor Torque')
grid on
subplot(3,1,2)
plot(simtime, xsim(:,11))
ylabel('$\Omega_2$ [rad/s]')
grid on
subplot(3,1,3)
plot(simtime, xsim(:,12))
xlabel('Time [s]')
ylabel('$\Omega_3$ [rad/s]')
grid on

% Body frame angular rate time history

% Define the Euler sequence matrix as a function of theta
E = @(th) [sin(th(3))*sin(th(2)) cos(th(3))  0;
           sin(th(2))*cos(th(3)) -sin(th(3)) 0;
           cos(th(2))            0           1];

% Map the Euler angle rates to the body angular velocity
w = zeros(length(simtime),3);
for k = 1:length(simtime)
    w(k,:) = (E(xsim(k,1:3))*xsim(k,7:9)')';
end

figure
subplot(3,1,1)
plot(simtime, w(:,1))
ylabel('$\omega_1$ [rad/s]')
title('Body Angular Rate Time History Zero Motor Torque')
grid on
subplot(3,1,2)
plot(simtime, w(:,2))
ylabel('$\omega_2$ [rad/s]')
grid on
subplot(3,1,3)
plot(simtime, w(:,3))
xlabel('Time [s]')
ylabel('$\omega_3$ [rad/s]')
grid on

% Alpha angular momentum time history

% Define palpha
pa = (w + xsim(:,10:12))*Jw;

% And plot
figure
subplot(3,1,1)
plot(simtime, pa(:,1))
ylabel('$p_{\alpha,1}$ [N m s]')
title('RW Generalized Momentum Time History Zero Motor Torque')
grid on
subplot(3,1,2)
plot(simtime, pa(:,2))
ylabel('$p_{\alpha,2}$ [N m s]')
grid on
subplot(3,1,3)
plot(simtime, pa(:,3))
xlabel('Time [s]')
ylabel('$p_{\alpha,3}$ [N m s]')
grid on

% Verification
C   = [300 -200 100];
dev = w + xsim(:,10:12) - C;          % deviation of the sum from C

figure
subplot(3,1,1)
plot(simtime, dev(:,1))
ylabel('$\omega_1+\Omega_1-C_1$')
title('Deviation of $\omega+\Omega$ from Initial Value')
grid on
subplot(3,1,2)
plot(simtime, dev(:,2))
ylabel('$\omega_2+\Omega_2-C_2$')
grid on
subplot(3,1,3)
plot(simtime, dev(:,3))
xlabel('Time [s]')
ylabel('$\omega_3+\Omega_3-C_3$')
grid on

% Export figures to ../Report
outdir = fullfile('..', 'Report');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
fig = figure(1);
exportgraphics(fig, fullfile(outdir, 'b2_theta.png'), 'Resolution', 300);
fig = figure(2);
exportgraphics(fig, fullfile(outdir, 'b2_wheelrate.png'), 'Resolution', 300);
fig = figure(3);
exportgraphics(fig, fullfile(outdir, 'b2_omega.png'), 'Resolution', 300);
fig = figure(4);
exportgraphics(fig, fullfile(outdir, 'b3_palpha.png'), 'Resolution', 300);
fig = figure(5);
exportgraphics(fig, fullfile(outdir, 'b4_check.png'), 'Resolution', 300);