% Recompute the initial conditions:
xkep = [6978; 0.01; deg2rad(45); 0; 0; 0];
mu = 398600;
x0 = kep2cart(xkep, mu);
% Compute the orbital period
a = 6978;
P = 2*pi*sqrt((a^3)/mu)*10;
% Run the simulation for P
tol = 10^(-12);
opt = odeset('RelTol', tol, 'AbsTol', tol);
[simtime, x] = ode113(@(t,x) twoBodyGrav(x ,mu), [0 P], x0, opt);

% Extract initial values
r = x0(1:3);
v = x0(4:6);
rnorm = norm(r);
vnorm = norm(v);

% Find specific energy
E = ((vnorm^2)/(2)) - ((mu)/(rnorm));

% Find specific angular momentum
h = cross(r,v);

% Find eccentricity vector
e = ((cross(v,h))/(mu))-(r/rnorm);

% Plot 1: trajectory only
figure
plot3(x(:,1), x(:,2), x(:,3))
axis equal
grid on
xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('J2 perturbed two-body trajectory (10 orbital periods)')

% Plot 2: trajectory with Earth
figure
plot3(x(:,1), x(:,2), x(:,3))
axis equal
hold on
grid on

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
% Apply texture to surface
surf(xs, ys, zs, 'FaceColor', 'texturemap', 'CData', flipud(earth_img), 'EdgeColor', 'none');

xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
title('J2 perturbed two-body trajectory with Earth (10 orbital periods)')

% Set LaTeX as default interpreter for all text
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLineLineWidth', 1.5);

% Compute time histories of E, h, e over the full simulation
N = size(x,1);
Et = zeros(N,1);
ht = zeros(N,3);
et = zeros(N,3);
dot_hr = zeros(N,1);
dot_hv = zeros(N,1);
dot_eh = zeros(N,1);

for k = 1:N
    rk = x(k,1:3)';
    vk = x(k,4:6)';
    rknorm = norm(rk);

    Et(k) = norm(vk)^2/2 - mu/rknorm;
    ht(k,:) = cross(rk, vk)';
    et(k,:) = (cross(vk, ht(k,:)')/mu - rk/rknorm)';

    dot_hr(k) = dot(ht(k,:)', rk);
    dot_hv(k) = dot(ht(k,:)', vk);
    dot_eh(k) = dot(et(k,:)', ht(k,:)');
end

% Plot 3: specific energy time history
figure
plot(simtime, Et)
xlabel('Time [s]')
ylabel('$\mathcal{E}$ [km$^2$/s$^2$]')
title('Specific energy time history')
grid on

% Plot 4: specific angular momentum time history (each component on its own subplot)
figure
subplot(3,1,1)
plot(simtime, ht(:,1))
ylabel('$h_x$ [km$^2$/s]')
title('Specific angular momentum time history')
grid on

subplot(3,1,2)
plot(simtime, ht(:,2))
ylabel('$h_y$ [km$^2$/s]')
grid on

subplot(3,1,3)
plot(simtime, ht(:,3))
xlabel('Time [s]')
ylabel('$h_z$ [km$^2$/s]')
grid on

% Plot 5: eccentricity vector time history (each component on its own subplot)
figure
subplot(3,1,1)
plot(simtime, et(:,1))
ylabel('$e_x$ [-]')
title('Eccentricity vector time history')
grid on

subplot(3,1,2)
plot(simtime, et(:,2))
ylabel('$e_y$ [-]')
grid on

subplot(3,1,3)
plot(simtime, et(:,3))
xlabel('Time [s]')
ylabel('$e_z$ [-]')
grid on

% Report orthogonality checks over the full simulation
fprintf('Max |h.r| over simulation: %e\n', max(abs(dot_hr)))
fprintf('Max |h.v| over simulation: %e\n', max(abs(dot_hv)))
fprintf('Max |e.h| over simulation: %e\n', max(abs(dot_eh)))

function xdot = twoBodyGrav(x, mu)
    J2 = 1.08263*10^(-3);
    r0 = 6378.137;

    r=x(1:3);
    v=x(4:6);
    rnorm = norm(r);
    aJ2 = ((-3*mu*J2*(r0^2))/(2*(rnorm^5))) * [ ...
    (1-(5*((r(3)^2)/(rnorm^2))))*r(1); ...
    (1-(5*((r(3)^2)/(rnorm^2))))*r(2); ...
    (3-(5*((r(3)^2)/(rnorm^2))))*r(3) ];
    xdot = [v; r*((-mu)/(rnorm^3))]+ [zeros(3,1); aJ2];
end