% Define the time variables and the beta m constant
t0 = 0;
tm = 60;
tf = 100;
t = linspace(t0,tf,10000);
T = t/tm;
bm = 60;

% Find the time histories of beta and its derivative
b = bm*(10*T.^3 - 15*T.^4 + 6*T.^5).*(t<=tm) + bm*(t>tm);
bd = bm/tm*(30*T.^2 - 60*T.^3 + 30*T.^4).*(t<=tm) + 0*(t>tm);

% Set LaTeX as default interpreter for all text
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLineLineWidth', 1.5);
set(groot, 'defaultAxesFontSize', 20);
set(groot, 'defaultTextFontSize', 20);
set(groot, 'defaultLegendFontSize', 20);

% Plot 1: beta and betadot time history (each quantity on its own subplot)
figure
subplot(2,1,1)
plot(t, b)
ylabel('$\beta$ [deg]')
title('Solar array angle time history')
grid on
subplot(2,1,2)
plot(t, bd)
xlabel('Time [s]')
ylabel('$\dot{\beta}$ [deg/s]')
grid on

% Export figures to ../Report
outdir = fullfile('..', 'Report');
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

% Beta and betadot from Plot 1
fig = figure(1);
exportgraphics(fig, fullfile(outdir, 'b1_beta.png'), 'Resolution', 300);