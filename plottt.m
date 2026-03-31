function plottt()
% =========================================================
% Square, bold, publication-quality figure
% PNG exported at 1200 dpi
% =========================================================

clc;
close all;

% -------------------- STYLE --------------------
AXIS_FS   = 28;
TICK_FS   = 24;
LEGEND_FS = 22;
LINE_W    = 4.5;
MARKER_SZ = 10;
AX_W      = 1.8;

set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

% -------------------- TIME --------------------
t = [0 4 8 12 16 24 32 40 48 56 64 72];

% -------------------- DATA --------------------
De0_01  = [1.0 0.045965866 0.034917319 0.023311509 0.013988117 0 0 0 0 0 0 0];
De0_025 = [1.0 0.046532867 0.03515301  0.023439505 0.014065955 0 0 0 0 0 0 0];
De0_05  = [1.0 0.048668343 0.036010081 0.023900508 0.01434501  5.62e-05 0 0 0 0 0 0];
De0_3   = [1.0 0.15308991  0.08872358  0.054278296 0.03301411  0.008364319 0 0 0 0 0 0];
De0_4   = [1.0 0.201239792 0.126211099 0.081472971 0.052351174 0.017994016 0 0 0 0 0 0];

% -------------------- FIGURE --------------------
fig = figure('Color','w','Position',[20 20 1000 1000]);
hold on;
box on;
grid on;

plot(t, De0_01,  '-o', 'LineWidth',LINE_W,'MarkerSize',MARKER_SZ, ...
    'MarkerFaceColor','auto','DisplayName','$L=0.1\,\mathrm{mm}$');
plot(t, De0_025, '-s', 'LineWidth',LINE_W,'MarkerSize',MARKER_SZ, ...
    'MarkerFaceColor','auto','DisplayName','$L=0.25\,\mathrm{mm}$');
plot(t, De0_05,  '-^', 'LineWidth',LINE_W,'MarkerSize',MARKER_SZ, ...
    'MarkerFaceColor','auto','DisplayName','$L=0.5\,\mathrm{mm}$');
plot(t, De0_3,   '-d', 'LineWidth',LINE_W,'MarkerSize',MARKER_SZ, ...
    'MarkerFaceColor','auto','DisplayName','$L=3\,\mathrm{mm}$');
plot(t, De0_4,   '-v', 'LineWidth',LINE_W,'MarkerSize',MARKER_SZ, ...
    'MarkerFaceColor','auto','DisplayName','$L=4\,\mathrm{mm}$');

% -------------------- AXES --------------------
xlabel('Time (hours)','FontSize',AXIS_FS,'Interpreter','latex');
ylabel('Total Sobol index of $D_{e0}$','FontSize',AXIS_FS,'Interpreter','latex');

set(gca, ...
    'FontSize',TICK_FS, ...
    'LineWidth',AX_W, ...
    'TickLength',[0.02 0.02], ...
    'TickLabelInterpreter','latex');

xlim([0 72]);
ylim([0 1.05]);

ax = gca;
ax.GridAlpha = 0.12;

% -------------------- LEGEND --------------------
lgd = legend('Location','northeast');
set(lgd, ...
    'FontSize',LEGEND_FS, ...
    'Box','off', ...
    'Interpreter','latex');

% -------------------- SAVE PNG ONLY --------------------
if ispc
    outdir = fullfile(getenv('USERPROFILE'),'Downloads');
else
    outdir = fullfile(getenv('HOME'),'Downloads');
end
if ~exist(outdir,'dir')
    mkdir(outdir);
end

pngfile = fullfile(outdir,'De0_square_1200dpi.png');

exportgraphics(fig, pngfile, 'Resolution', 1200, 'BackgroundColor', 'white');

fprintf('Saved PNG: %s\n', pngfile);

end