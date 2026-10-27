% Define the beta angle function
bm = deg2rad(60);
tm = 60;
bt = @(t) bm*(10*(t/tm).^3-15*(t/tm).^4+6*(t/tm).^5).*(t<=tm)+bm*(t>tm);

% Define the time vector
t = linspace(0,60,10000);

% Define the torque magnitude
T = 1;

% Compute the panel normal torque at each time
u = T*[sin(bt(t)); zeros(1,length(t)); cos(bt(t))];


% Set LaTeX as default interpreter for all text
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultAxesFontSize', 20);
set(groot, 'defaultTextFontSize', 20);
set(groot, 'defaultLegendFontSize', 20);

% RW motor torque time history
figure
subplot(3,1,1)
plot(t, u(1,:))
ylabel('$u_1$ [N m]')
title('RW Motor Torque Time History')
grid on
subplot(3,1,2)
plot(t, u(2,:))
ylabel('$u_2$ [N m]')
grid on
subplot(3,1,3)
plot(t, u(3,:))
xlabel('Time [s]')
ylabel('$u_3$ [N m]')
grid on

% Export figures to ../Report
outdir = fullfile('..', 'Report');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
fig = figure(1);
exportgraphics(fig, fullfile(outdir, 'c1_torque.png'), 'Resolution', 300);