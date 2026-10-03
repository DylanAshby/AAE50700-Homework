% Set LaTeX as default interpreter for all text
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultAxesFontSize', 20);
set(groot, 'defaultTextFontSize', 20);
set(groot, 'defaultLegendFontSize', 20);

% Define the relevant values
a = 6628;
e = 0.01;
mu = 398600;
J2 = 1.08263*10^(-3);
ro = 6378.137;

% Define the drag parameters, converted to km-based units
CD = 2.2;
A = 0.01*(1/1000)^2;            % A/m [km^2/kg]
rho0 = 1.454*10^(-13)*(1000)^3; % [kg/km^3]
H = 60;                         % [km]

% Compute the propagation time of 5 days in seconds
tf = 5*24*60*60;

% Compute r and v magnitude at periapsis
p = a*(1-e^2);
r = p/(1+(e*cos(0)));
v = sqrt(((2*mu)/r)-(mu/a));

% Define the vector forms in RTN
rB = [r; 0; 0];
vB = [0; v; 0];

% Define the DCM
NB = [1 0 0; 0 cosd(45) -sind(45); 0 sind(45) cosd(45)];

% Then move rB and vB to ECI frame
rN = NB*rB;
vN = NB*vB;

% Compute spherical initial conditions
x = rN(1); y = rN(2); z = rN(3); xd = vN(1); yd = vN(2); zd = vN(3);
rS = [r; atan2d(y,x); atan2d(z,sqrt(x^2+y^2))];
vS = [((x*xd)+(y*yd)+(z*zd))/r;
      rad2deg(((x*yd)-(y*xd))/((x^2)+(y^2)));
      rad2deg(((((x^2)+(y^2))*zd)-(z*((x*xd)+(y*yd))))/(sqrt((x^2)+(y^2))*(r^2)))];

% Build the integrator initial state in radians
x0 = [rS(1); deg2rad(rS(2)); deg2rad(rS(3)); vS(1); deg2rad(vS(2)); deg2rad(vS(3))];

% Run the simulation for both models on a shared time vector
tspan = 0:10:tf;
tol = 10^(-12);
opt = odeset('RelTol', tol, 'AbsTol', tol);
[simtime, xnodrag] = ode45(@(t,x) J2nodrag(x, mu, J2, ro), tspan, x0, opt);
[~, xdrag] = ode45(@(t,x) J2drag(x, mu, J2, ro, CD, A, rho0, H), tspan, x0, opt);
[~, xaug] = ode45(@(t,x) J2dragaug(x, mu, J2, ro, CD, A, rho0, H), tspan, [x0; 0], opt);

% Convert spherical states to Cartesian positions
rnodrag = [xnodrag(:,1).*cos(xnodrag(:,3)).*cos(xnodrag(:,2)), ...
           xnodrag(:,1).*cos(xnodrag(:,3)).*sin(xnodrag(:,2)), ...
           xnodrag(:,1).*sin(xnodrag(:,3))];
rdrag = [xdrag(:,1).*cos(xdrag(:,3)).*cos(xdrag(:,2)), ...
         xdrag(:,1).*cos(xdrag(:,3)).*sin(xdrag(:,2)), ...
         xdrag(:,1).*sin(xdrag(:,3))];

% Load Earth texture
earth_img = imread('2k_earth_daymap.jpg');

% Generate sphere geometry
earth_radius = 6378.137;
earth_center = [0; 0; 0];
n = 3000;                       % resolution
[xs, ys, zs] = sphere(n);
xs = xs * earth_radius + earth_center(1);
ys = ys * earth_radius + earth_center(2);
zs = zs * earth_radius + earth_center(3);

% Compute the specific Jacobi integral and z angular momentum for both models
Jnodrag = (xnodrag(:,4).^2)/2 + ((xnodrag(:,1).^2).*(cos(xnodrag(:,3)).^2).*(xnodrag(:,5).^2))/2 + ...
          ((xnodrag(:,1).^2).*(xnodrag(:,6).^2))/2 - mu./xnodrag(:,1) + ...
          ((3*mu*J2*(ro^2))./(2*(xnodrag(:,1).^3))).*(sin(xnodrag(:,3)).^2) - ...
          (mu*J2*(ro^2))./(2*(xnodrag(:,1).^3));
Jdrag = (xdrag(:,4).^2)/2 + ((xdrag(:,1).^2).*(cos(xdrag(:,3)).^2).*(xdrag(:,5).^2))/2 + ...
        ((xdrag(:,1).^2).*(xdrag(:,6).^2))/2 - mu./xdrag(:,1) + ...
        ((3*mu*J2*(ro^2))./(2*(xdrag(:,1).^3))).*(sin(xdrag(:,3)).^2) - ...
        (mu*J2*(ro^2))./(2*(xdrag(:,1).^3));
Hznodrag = (xnodrag(:,1).^2).*xnodrag(:,5).*(cos(xnodrag(:,3)).^2);
Hzdrag = (xdrag(:,1).^2).*xdrag(:,5).*(cos(xdrag(:,3)).^2);

% Compute the normalized changes from the initial values
dJnodrag = (Jnodrag - Jnodrag(1))/abs(Jnodrag(1));
dJdrag = (Jdrag - Jdrag(1))/abs(Jdrag(1));
dHznodrag = (Hznodrag - Hznodrag(1))/abs(Hznodrag(1));
dHzdrag = (Hzdrag - Hzdrag(1))/abs(Hzdrag(1));

% Compute the specific Jacobi integral for the augmented run
Jaug = (xaug(:,4).^2)/2 + ((xaug(:,1).^2).*(cos(xaug(:,3)).^2).*(xaug(:,5).^2))/2 + ...
       ((xaug(:,1).^2).*(xaug(:,6).^2))/2 - mu./xaug(:,1) + ...
       ((3*mu*J2*(ro^2))./(2*(xaug(:,1).^3))).*(sin(xaug(:,3)).^2) - ...
       (mu*J2*(ro^2))./(2*(xaug(:,1).^3));

% Compute the normalized work-balance residual
resid = (Jaug - Jaug(1) - xaug(:,7))/abs(Jaug(1));

% Plot 1: J2 only trajectory
figure
plot3(rnodrag(:,1), rnodrag(:,2), rnodrag(:,3))
axis equal
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('$J_2$ only trajectory (5 days)')

% Plot 2: J2 only trajectory with Earth
figure
plot3(rnodrag(:,1), rnodrag(:,2), rnodrag(:,3))
axis equal
hold on
grid on
surf(xs, ys, zs, 'FaceColor', 'texturemap', 'CData', flipud(earth_img), 'EdgeColor', 'none');
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('$J_2$ only trajectory with Earth (5 days)')

% Plot 3: J2 with drag trajectory
figure
plot3(rdrag(:,1), rdrag(:,2), rdrag(:,3))
axis equal
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('$J_2$ with drag trajectory (5 days)')

% Plot 4: J2 with drag trajectory with Earth
figure
plot3(rdrag(:,1), rdrag(:,2), rdrag(:,3))
axis equal
hold on
grid on
surf(xs, ys, zs, 'FaceColor', 'texturemap', 'CData', flipud(earth_img), 'EdgeColor', 'none');
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('$J_2$ with drag trajectory with Earth (5 days)')

% Plot 5: position difference time history (each component on its own subplot)
dr = rdrag - rnodrag;
tdays = simtime/(24*60*60);
figure
subplot(3,1,1)
plot(tdays, dr(:,1))
ylabel('$\Delta x$ [km]')
title('Position difference time history (drag minus $J_2$ only)')
grid on
subplot(3,1,2)
plot(tdays, dr(:,2))
ylabel('$\Delta y$ [km]')
grid on
subplot(3,1,3)
plot(tdays, dr(:,3))
xlabel('Time [days]')
ylabel('$\Delta z$ [km]')
grid on

% Plot 6: J2 only normalized changes (each quantity on its own subplot)
figure
subplot(2,1,1)
plot(tdays, dJnodrag)
ylabel('$\frac{\mathcal{J}(t)-\mathcal{J}(0)}{|\mathcal{J}(0)|}$')
title('$J_2$ only normalized changes')
grid on
subplot(2,1,2)
plot(tdays, dHznodrag)
xlabel('Time [days]')
ylabel('$\frac{H_z(t)-H_z(0)}{|H_z(0)|}$')
grid on

% Plot 7: J2 with drag normalized changes (each quantity on its own subplot)
figure
subplot(2,1,1)
plot(tdays, dJdrag)
ylabel('$\frac{\mathcal{J}(t)-\mathcal{J}(0)}{|\mathcal{J}(0)|}$')
title('$J_2$ with drag normalized changes')
grid on
subplot(2,1,2)
plot(tdays, dHzdrag)
xlabel('Time [days]')
ylabel('$\frac{H_z(t)-H_z(0)}{|H_z(0)|}$')
grid on

% Plot 8: normalized work-balance residual
figure
plot(tdays, resid)
xlabel('Time [days]')
ylabel('$\frac{\mathcal{J}(t)-\mathcal{J}(0)-W_D(t)}{|\mathcal{J}(0)|}$')
title('Normalized work-balance residual')
grid on

% Export figures to ../Report
outdir = fullfile('..', 'Report');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% J2 only: 3D, edge-on (X-Z), and side-on (X-Y) views from Plot 1
fig = figure(1); ax = fig.CurrentAxes;
view(ax, 3);
exportgraphics(fig, fullfile(outdir, 'c2_J2only_3d.png'), 'Resolution', 300);
view(ax, -90, 0);
exportgraphics(fig, fullfile(outdir, 'c2_J2only_edgeon.png'), 'Resolution', 300);
view(ax, 0, 90);
exportgraphics(fig, fullfile(outdir, 'c2_J2only_sideon.png'), 'Resolution', 300);
view(ax, 3);

% J2 only with Earth from Plot 2
fig = figure(2);
exportgraphics(fig, fullfile(outdir, 'c2_J2only_earth.png'), 'Resolution', 300);

% J2 with drag: 3D, edge-on (X-Z), and side-on (X-Y) views from Plot 3
fig = figure(3); ax = fig.CurrentAxes;
view(ax, 3);
exportgraphics(fig, fullfile(outdir, 'c2_J2drag_3d.png'), 'Resolution', 300);
view(ax, -90, 0);
exportgraphics(fig, fullfile(outdir, 'c2_J2drag_edgeon.png'), 'Resolution', 300);
view(ax, 0, 90);
exportgraphics(fig, fullfile(outdir, 'c2_J2drag_sideon.png'), 'Resolution', 300);
view(ax, 3);

% J2 with drag with Earth from Plot 4
fig = figure(4);
exportgraphics(fig, fullfile(outdir, 'c2_J2drag_earth.png'), 'Resolution', 300);

% Position difference from Plot 5
fig = figure(5);
exportgraphics(fig, fullfile(outdir, 'c2_difference.png'), 'Resolution', 300);

% J2 only normalized changes from Plot 6
fig = figure(6);
exportgraphics(fig, fullfile(outdir, 'c3_J2only_normalized.png'), 'Resolution', 300);

% J2 with drag normalized changes from Plot 7
fig = figure(7);
exportgraphics(fig, fullfile(outdir, 'c3_J2drag_normalized.png'), 'Resolution', 300);

% Work-balance residual from Plot 8
fig = figure(8);
exportgraphics(fig, fullfile(outdir, 'c4_residual.png'), 'Resolution', 300);