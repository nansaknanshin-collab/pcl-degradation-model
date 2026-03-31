function DimensionalPlot()
% =========================================================
% DIMENSIONAL JOINT FIT (mm + hours)
% Save this file as: DimensionalPlot.m
% Run using: DimensionalPlot
% =========================================================
rng(1);

% -------------------- GLOBAL PLOT SETTINGS --------------------
FS = 20;           % font size for axis labels, tick labels, and legends
EXPORT_DPI = 1200; % export resolution

% -------------------- SETTINGS --------------------
MIN_ERR_PCT = 0.1;   % error floor (%)
Ebulk       = 1;     % Dirichlet boundary enzyme at x=0
Nx          = 25;    % grid points
L_mm        = 0.25;  % half-thickness domain (mm)
dx          = L_mm/(Nx-1);
alphaE      = 1.0;

% -------------------- Parameter bounds --------------------
RATE_LB = 1e-6;
RATE_UB = 50;
De0_LB  = 1e-12;
De0_UB  = 50;

% -------------------- MCMC settings --------------------
nit       = 500;
nthin     = 5;
nchain    = 5;
n_walkers = nchain * 8;
gamma_val = [];
jitter    = [];

% -------------------- Fixed initial polymer fractions --------------------
C0 = 0.302;
A0 = 0.698;

% =========================================================
% DATASET 1: Weight loss (hours, %, sigma%)
% =========================================================
P_raw = [ ...
    72, 87.56, 0.1; ...
    64, 86.12, 3.0; ...
    56, 86.06, 1.5; ...
    48, 86.03, 1.4; ...
    40, 84.80, 1.4; ...
    32, 83.75, 1.5; ...
    24, 81.47, 1.4; ...
    20, 80.47, 3.0; ...
    16, 75.67, 1.5; ...
    12, 66.66, 1.2; ...
     8, 54.85, 1.5; ...
     4, 38.12, 1.2; ...
     0,  0.00, 0.0  ...
];
[timeP, ordP] = sort(P_raw(:,1));
dataP = P_raw(ordP,2);
errsP = max(P_raw(ordP,3), MIN_ERR_PCT);

% =========================================================
% DATASET 2: DSC crystallinity (hours, %, sigma%)
% =========================================================
shi_lipase = [ ...
     0, 30.2, 0.9; ...
     4, 21.6, 0.8; ...
     8, 17.1, 0.6; ...
    16,  7.2, 0.9; ...
    24,  6.9, 0.6; ...
    40,  0.0, 0.0; ...
    64,  0.0, 0.0  ...
];
[timeC, ordC] = sort(shi_lipase(:,1));
dataXC = shi_lipase(ordC,2);
errsXC = max(shi_lipase(ordC,3), MIN_ERR_PCT);

% =========================================================
% INITIAL GUESS
% =========================================================
k1_init    = 0.651327;
km1_init   = 0.474921;
k3_init    = 0.300349;
km3_init   = 0.0282859;
kconv_init = 0.575375;
kdegC_init = 0.00116811;
kdegA_init = 0.00492959;
De0_init   = 1.65766;

pml = [k1_init, km1_init, k3_init, km3_init, kconv_init, ...
       kdegC_init, kdegA_init, De0_init];

fprintf('\n=== SHI DATA: JOINT FIT — DIMENSIONAL (mm + hours) ===\n');
fprintf('L = %.6g mm, Nx = %d, dx = %.6g mm\n', L_mm, Nx, dx);
fprintf('Ebulk = %.3g, alphaE = %.3g\n', Ebulk, alphaE);
fprintf('Bounds: rates [%.2e, %.2g], De0 [%.2e, %.2g]\n', ...
    RATE_LB, RATE_UB, De0_LB, De0_UB);
fprintf('NOTE: Constraint enforced: kdegA > kdegC.\n\n');

% -------------------- Run MCMC (DE-MCz) --------------------
[samples_log, ~] = run_block_joint( ...
    timeP, dataP, errsP, ...
    timeC, dataXC, errsXC, ...
    C0, A0, Ebulk, Nx, dx, alphaE, ...
    pml, nit, nthin, n_walkers, gamma_val, jitter, ...
    RATE_LB, RATE_UB, De0_LB, De0_UB);

if isempty(samples_log)
    error('No posterior samples saved. Likelihood may be -inf everywhere.');
end

samples = exp(samples_log);
samples = samples(all(isfinite(samples),2),:);

if isempty(samples)
    error('All posterior samples became non-finite after conversion.');
end

best_params = mean(samples, 1);

param_labels = {'$k_1$','$k_{-1}$','$k_3$','$k_{-3}$', ...
                '$k_{\mathrm{conv}}$','$k_{\mathrm{deg},C}$', ...
                '$k_{\mathrm{deg},A}$','$D_{e0}$'};

fprintf('Posterior mean parameters (physical):\n');
for i = 1:numel(best_params)
    if i == 8
        fprintf('%s = %.6g  [mm^2/h]\n', param_labels{i}, best_params(i));
    else
        fprintf('%s = %.6g  [1/h]\n', param_labels{i}, best_params(i));
    end
end
fprintf('\nSaved posterior samples (finite): %d\n', size(samples,1));

% -------------------- CI plots --------------------
figure('Color','w','Position',[100 100 1000 500]);
plot_posterior_CI_main_rates(samples, param_labels, FS);
save_fig('Figure 9.png', EXPORT_DPI);

figure('Color','w','Position',[100 100 600 500]);
plot_posterior_CI_small_rates(samples, param_labels, FS);
save_fig('Figure 10.png', EXPORT_DPI);

figure('Color','w','Position',[100 100 400 500]);
plot_posterior_CI_De0(samples, param_labels, FS);
save_fig('Figure 11.png', EXPORT_DPI);

% -------------------- Prediction uncertainty --------------------
Q95_WL = compute_Q95_WL(samples, best_params, timeP, C0, A0, Ebulk, Nx, dx, alphaE);
Q95_XC = compute_Q95_XC(samples, best_params, timeC, C0, A0, Ebulk, Nx, dx, alphaE);
fprintf('\nQ_{0.95} (WL) = %.6f\n', Q95_WL);
fprintf('Q_{0.95} (Xc) = %.6f\n', Q95_XC);

tmin  = min([timeP(:); timeC(:)]);
tmax  = max([timeP(:); timeC(:)]);
tgrid = linspace(tmin, tmax, 600);

[predWL, predXC] = compute_predictions_joint(samples, C0, A0, Ebulk, Nx, dx, alphaE, tgrid);

ciWL_low  = prctile(predWL,  2.5, 1);
ciWL_high = prctile(predWL, 97.5, 1);
ciXC_low  = prctile(predXC,  2.5, 1);
ciXC_high = prctile(predXC, 97.5, 1);
meanWL = mean(predWL, 1);
meanXC = mean(predXC, 1);

figure('Color','w','Position',[100 100 900 650]); hold on; grid on;
xlabel('Time (hours)','FontSize',FS);
ylabel('Weight Loss (%)','FontSize',FS);
errorbar(timeP, dataP, errsP, errsP, 'o', ...
    'MarkerSize',6, 'Color','k', 'MarkerFaceColor','k', ...
    'LineStyle','none', 'DisplayName','PCL-Lipase data');
fill([tgrid fliplr(tgrid)], [ciWL_low fliplr(ciWL_high)], ...
    [0 0.7 0], 'EdgeColor','none', 'FaceAlpha',0.35, ...
    'DisplayName','Prediction uncertainty');
plot(tgrid, meanWL, 'r-', 'LineWidth',3, 'DisplayName','Model fit');
set(gca,'FontSize',FS);
legend('Location','southeast','FontSize',FS,'Box','off');
save_fig('Shi_WL_PredictionUncertainty.png', EXPORT_DPI);

figure('Color','w','Position',[100 100 900 650]); hold on; grid on;
xlabel('Time (hours)','FontSize',FS);
ylabel('$\chi_c$ (\%)','FontSize',FS,'Interpreter','latex');
errorbar(timeC, dataXC, errsXC, errsXC, 'o', ...
    'MarkerSize',6, 'Color','k', 'MarkerFaceColor','k', ...
    'LineStyle','none', 'DisplayName','$\chi_c$ data');
fill([tgrid fliplr(tgrid)], [ciXC_low fliplr(ciXC_high)], ...
    [0 0.7 0], 'EdgeColor','none', 'FaceAlpha',0.35, ...
    'DisplayName','Prediction uncertainty');
plot(tgrid, meanXC, 'r-', 'LineWidth',3, 'DisplayName','Model fit');
set(gca,'FontSize',FS);
lgd = legend('Location','northeast','FontSize',FS,'Box','off');
set(lgd,'Interpreter','latex');
save_fig('Shi_Xc_PredictionUncertainty.png', EXPORT_DPI);

% -------------------- Pairwise posterior --------------------
plot_pairwise(samples, best_params, param_labels, FS);
save_fig('Shi_PairwisePosterior.png', EXPORT_DPI);

fprintf('\nDone.\n');
end

function [saved_states, saved_logprob] = run_block_joint( ...
    timeP, dataP, errsP, ...
    timeC, dataXC, errsXC, ...
    C0, A0, Ebulk, Nx, dx, alphaE, ...
    pml, nit, nthin, n_walkers, gamma_val, jitter, ...
    RATE_LB, RATE_UB, De0_LB, De0_UB)

lower = [RATE_LB, RATE_LB, RATE_LB, RATE_LB, RATE_LB, RATE_LB, RATE_LB, De0_LB];
upper = [RATE_UB, RATE_UB, RATE_UB, RATE_UB, RATE_UB, RATE_UB, RATE_UB, De0_UB];
lb = log(lower);
ub = log(upper);

d = numel(pml);

logpost = @(theta) logposterior_joint(theta, lb, ub, ...
    timeP, dataP, errsP, timeC, dataXC, errsXC, C0, A0, Ebulk, Nx, dx, alphaE);

init_pos = nan(n_walkers, d);
for w = 1:n_walkers
    ok = false;
    for tries = 1:200
        cand = log(pml) + 0.1*randn(1,d);
        cand = min(max(cand, lb), ub);
        if isfinite(logpost(cand))
            init_pos(w,:) = cand;
            ok = true;
            break;
        end
    end
    if ~ok
        init_pos(w,:) = min(max(log(pml), lb), ub);
    end
end

Z = init_pos;
current_logprob = [];

[states_out, saved_logprob] = DEMCz( ...
    logpost, init_pos, Z, nit, nthin, {1:d}, current_logprob, [], gamma_val, jitter);

mask = isfinite(saved_logprob) & all(states_out >= lb & states_out <= ub, 2);
saved_states  = states_out(mask,:);
saved_logprob = saved_logprob(mask,:);
end

function logp = logposterior_joint(theta, lb, ub, ...
    timeP, dataP, errsP, timeC, dataXC, errsXC, ...
    C0, A0, Ebulk, Nx, dx, alphaE)

if any(theta < lb) || any(theta > ub)
    logp = -inf;
    return;
end

phys = exp(theta);

kdegC = phys(6);
kdegA = phys(7);
if kdegA <= kdegC
    logp = -inf;
    return;
end

try
    t_all = sort(unique([timeP(:); timeC(:)]));
    [WL_all, XC_all] = simulate_PDE_observables(t_all, phys, C0, A0, Ebulk, Nx, dx, alphaE);

    WL_model = interp1(t_all, WL_all, timeP(:), 'linear', 'extrap');
    XC_model = interp1(t_all, XC_all, timeC(:), 'linear', 'extrap');

    if any(~isfinite(WL_model)) || any(~isfinite(XC_model))
        logp = -inf;
        return;
    end

    rWL = (WL_model - dataP(:)) ./ errsP(:);
    rXC = (XC_model - dataXC(:)) ./ errsXC(:);

    logp = -0.5*(sum(rWL.^2) + sum(rXC.^2)) ...
         - (sum(log(errsP(:))) + sum(log(errsXC(:))));

    if ~isfinite(logp)
        logp = -inf;
    end
catch
    logp = -inf;
end
end

function [WL, XC] = simulate_PDE_observables(t_eval, params, C0, A0, Ebulk, Nx, dx, alphaE)
k1    = params(1);
km1   = params(2);
k3    = params(3);
km3   = params(4);
kconv = params(5);
kdegC = params(6);
kdegA = params(7);
De0   = params(8);

E0 = zeros(Nx,1);
E0(1) = Ebulk;
C  = C0*ones(Nx,1);
A  = A0*ones(Nx,1);
EC = zeros(Nx,1);
EA = zeros(Nx,1);
P  = zeros(Nx,1);
y0 = [E0; C; A; EC; EA; P];

t_eval = sort(t_eval(:));
opts = odeset('RelTol',1e-6,'AbsTol',1e-9,'MaxStep',0.5);
[~, Y] = ode15s(@rhs_mol, t_eval, y0, opts);

if any(~isfinite(Y(:)))
    error('Non-finite ODE solution');
end

Cmat = Y(:, (1*Nx+1):(2*Nx));
Amat = Y(:, (2*Nx+1):(3*Nx));

L_mm_local = (Nx-1)*dx;
x = linspace(0, L_mm_local, Nx);

M0 = trapz(x, Cmat(1,:) + Amat(1,:)) / L_mm_local;
M  = trapz(x, Cmat + Amat, 2) / L_mm_local;
WL = 100*(1 - M./max(M0,1e-12));

denom = trapz(x, Cmat + Amat, 2) / L_mm_local;
Cbar  = trapz(x, Cmat, 2) / L_mm_local;
XC = 100*(Cbar ./ max(denom,1e-12));

    function dydt = rhs_mol(~, y)
        E  = y(1:Nx);
        C  = y((1*Nx+1):(2*Nx));
        A  = y((2*Nx+1):(3*Nx));
        EC = y((3*Nx+1):(4*Nx));
        EA = y((4*Nx+1):(5*Nx));
        P  = y((5*Nx+1):(6*Nx));

        E(1) = Ebulk;

        denom0_local = max(C0 + A0, 1e-12);
        phi = 1 - (C + A) ./ denom0_local;
        phi = min(max(phi, -0.5), 2.0);

        De = De0 .* (1 + alphaE .* phi);
        De = max(De, 1e-30);

        diffE = diffusion_dirichlet_neumann_varD(E, De, Ebulk, dx);

        dE  = diffE + (-(k1.*E.*C) + (km1.*EC) ...
                       -(k3.*E.*A) + (km3.*EA) ...
                       + (kconv.*EC) + (kdegC.*EC) + (kdegA.*EA));

        dC  = -(k1.*E.*C) + (km1.*EC);
        dA  = -(k3.*E.*A) + (km3.*EA) + (kconv.*EC);
        dEC =  (k1.*E.*C) - (km1.*EC) - (kconv.*EC) - (kdegC.*EC);
        dEA =  (k3.*E.*A) - (km3.*EA) - (kdegA.*EA);
        dP  =  (kdegC.*EC) + (kdegA.*EA);

        dE(1) = 0;
        dydt = [dE; dC; dA; dEC; dEA; dP];
    end
end

function div = diffusion_dirichlet_neumann_varD(u, D, u_left, dx)
u = u(:);
D = D(:);
Nx_local = numel(u);
u(1) = u_left;

flux = zeros(Nx_local+1,1);
D_iface = 0.5*(D(1:Nx_local-1) + D(2:Nx_local));
flux(2:Nx_local) = D_iface .* (u(2:Nx_local) - u(1:Nx_local-1)) / dx;
flux(Nx_local+1) = 0;

div = zeros(Nx_local,1);
div(1) = 0;
div(2:Nx_local-1) = (flux(3:Nx_local) - flux(2:Nx_local-1)) / dx;
div(Nx_local) = 2*(flux(Nx_local+1) - flux(Nx_local)) / dx;
end

function [predWL, predXC] = compute_predictions_joint(samples, C0, A0, Ebulk, Nx, dx, alphaE, tgrid)
nsamp = size(samples,1);
npnts = numel(tgrid);
predWL = zeros(nsamp, npnts);
predXC = zeros(nsamp, npnts);

for i = 1:nsamp
    [WL, XC] = simulate_PDE_observables(tgrid(:), samples(i,:), C0, A0, Ebulk, Nx, dx, alphaE);
    predWL(i,:) = WL(:).';
    predXC(i,:) = XC(:).';
end
end

function Q95 = compute_Q95_WL(samples, best_params, time, C0, A0, Ebulk, Nx, dx, alphaE)
nsamp = size(samples,1);
[WL_best, ~] = simulate_PDE_observables(time(:), best_params, C0, A0, Ebulk, Nx, dx, alphaE);
Q = nan(nsamp,1);

for i = 1:nsamp
    [WL_i, ~] = simulate_PDE_observables(time(:), samples(i,:), C0, A0, Ebulk, Nx, dx, alphaE);
    valid = (WL_i > 1e-8) & (WL_best > 1e-8);
    if any(valid)
        ld = log2(WL_i(valid) ./ WL_best(valid));
        Q(i) = sum(ld.^2);
    end
end

Q = Q(isfinite(Q));
Q95 = prctile(Q, 95);
end

function Q95 = compute_Q95_XC(samples, best_params, time, C0, A0, Ebulk, Nx, dx, alphaE)
nsamp = size(samples,1);
[~, XC_best] = simulate_PDE_observables(time(:), best_params, C0, A0, Ebulk, Nx, dx, alphaE);
Q = nan(nsamp,1);

for i = 1:nsamp
    [~, XC_i] = simulate_PDE_observables(time(:), samples(i,:), C0, A0, Ebulk, Nx, dx, alphaE);
    valid = (XC_i > 1e-8) & (XC_best > 1e-8);
    if any(valid)
        ld = log2(XC_i(valid) ./ XC_best(valid));
        Q(i) = sum(ld.^2);
    end
end

Q = Q(isfinite(Q));
Q95 = prctile(Q, 95);
end

function plot_pairwise(samples, best_params, param_labels, FS)
n_dim = size(samples,2);
figure('Name','Pairwise Posterior','Color','w','Position',[100,100,1200,1000]);

n_rows = ceil(sqrt(n_dim-1));
n_cols = ceil((n_dim-1)/n_rows);

for i = 1:(n_dim-1)
    subplot(n_rows,n_cols,i); hold on; box on;
    scatter(samples(:,i), samples(:,i+1), 25, 'b', 'filled', ...
        'MarkerFaceAlpha',0.35);
    plot(best_params(i), best_params(i+1), 'r*', ...
        'MarkerSize',10,'LineWidth',1.5);
    xlabel(param_labels{i},   'Interpreter','latex','FontSize',FS);
    ylabel(param_labels{i+1}, 'Interpreter','latex','FontSize',FS);
    set(gca,'FontSize',FS);
    grid on;
end
end

function plot_posterior_CI_main_rates(samples_phys, param_names, FS)
samples_phys = samples_phys(all(isfinite(samples_phys),2),:);
means = mean(samples_phys,1);
lower = prctile(samples_phys,2.5,1);
upper = prctile(samples_phys,97.5,1);

rate_idx = [1 2 3 5];
cmap = lines(8);

cla; hold on; box on; grid on;
h_leg = gobjects(1,numel(rate_idx));

for ii = 1:numel(rate_idx)
    i = rate_idx(ii);
    rectangle('Position',[ii-0.18, lower(i), 0.36, max(upper(i)-lower(i),0)], ...
        'FaceColor',cmap(i,:), 'EdgeColor',cmap(i,:), 'LineWidth',1.2);
    plot([ii-0.18 ii+0.18], [means(i) means(i)], 'r-', 'LineWidth',2);
    h_leg(ii) = patch(NaN, NaN, cmap(i,:), 'EdgeColor', cmap(i,:));
end

h_mean = plot(NaN, NaN, 'r-', 'LineWidth',2);

xlim([0 numel(rate_idx)+1]);
ylim([0 1.10*max(upper(rate_idx))]);

set(gca, 'XTick',1:numel(rate_idx), ...
    'XTickLabel',param_labels_subset(param_names, rate_idx), ...
    'TickLabelInterpreter','latex', 'FontSize',FS);

xlabel('Rate parameters','FontSize',FS,'Interpreter','latex');
ylabel('Parameter values ($\mathrm{h}^{-1}$)','FontSize',FS,'Interpreter','latex');

legend_labels = cell(1,numel(rate_idx));
for ii = 1:numel(rate_idx)
    legend_labels{ii} = sprintf('%s (95\\%% CI)', param_names{rate_idx(ii)});
end

legend([h_leg h_mean], [legend_labels {'Mean'}], ...
    'Location','northwest', 'Interpreter','latex', ...
    'FontSize',FS, 'NumColumns',1, 'Box','off');
end

function plot_posterior_CI_small_rates(samples_phys, param_names, FS)
samples_phys = samples_phys(all(isfinite(samples_phys),2),:);
means = mean(samples_phys,1);
lower = prctile(samples_phys,2.5,1);
upper = prctile(samples_phys,97.5,1);

small_idx = [4 6 7];
cmap = lines(8);

cla; hold on; box on; grid on;
h_leg = gobjects(1,numel(small_idx));

for ii = 1:numel(small_idx)
    i = small_idx(ii);
    rectangle('Position',[ii-0.18, lower(i), 0.36, max(upper(i)-lower(i),0)], ...
        'FaceColor',cmap(i,:), 'EdgeColor',cmap(i,:), 'LineWidth',1.2);
    plot([ii-0.18 ii+0.18], [means(i) means(i)], 'r-', 'LineWidth',2);
    h_leg(ii) = patch(NaN, NaN, cmap(i,:), 'EdgeColor', cmap(i,:));
end

h_mean = plot(NaN, NaN, 'r-', 'LineWidth',2);

xlim([0.5 numel(small_idx)+0.5]);
ylim([0 1.10*max(upper(small_idx))]);

set(gca, 'XTick',1:numel(small_idx), ...
    'XTickLabel',param_labels_subset(param_names, small_idx), ...
    'TickLabelInterpreter','latex', 'FontSize',FS);

xlabel('Rate parameters','FontSize',FS,'Interpreter','latex');
ylabel('Parameter values ($\mathrm{h}^{-1}$)','FontSize',FS,'Interpreter','latex');

legend_labels = cell(1,numel(small_idx));
for ii = 1:numel(small_idx)
    legend_labels{ii} = sprintf('%s (95\\%% CI)', param_names{small_idx(ii)});
end

legend([h_leg h_mean], [legend_labels {'Mean'}], ...
    'Location','northeast', 'Interpreter','latex', ...
    'FontSize',FS, 'NumColumns',1, 'Box','off');
end

function plot_posterior_CI_De0(samples_phys, param_names, FS)
samples_phys = samples_phys(all(isfinite(samples_phys),2),:);
means = mean(samples_phys,1);
lower = prctile(samples_phys,2.5,1);
upper = prctile(samples_phys,97.5,1);

De_idx = 8;
cmap = lines(8);

cla; hold on; box on; grid on;

rectangle('Position',[1-0.18, lower(De_idx), 0.36, max(upper(De_idx)-lower(De_idx),0)], ...
    'FaceColor',cmap(De_idx,:), 'EdgeColor',cmap(De_idx,:), 'LineWidth',1.2);

plot([1-0.18 1+0.18], [means(De_idx) means(De_idx)], 'r-', 'LineWidth',2);

h_ci   = patch(NaN, NaN, cmap(De_idx,:), 'EdgeColor', cmap(De_idx,:));
h_mean = plot(NaN, NaN, 'r-', 'LineWidth',2);

xlim([0.0 2.0]);
ylim([0 1.10*upper(De_idx)]);

set(gca, 'XTick',1, ...
    'XTickLabel',param_names(De_idx), ...
    'TickLabelInterpreter','latex', 'FontSize',FS);

xlabel('Diffusion parameter','FontSize',FS,'Interpreter','latex');
ylabel('Parameter values (mm/h)','FontSize',FS,'Interpreter','latex');

lgd = legend([h_ci h_mean], {'95\% CI','Mean'}, ...
    'Location','northeastoutside', 'Interpreter','latex', ...
    'FontSize',FS, 'NumColumns',1, 'Box','off');
set(lgd,'Color','none');
end

function labels = param_labels_subset(param_names, idx)
labels = param_names(idx);
end

function save_fig(fname, dpi)
if nargin < 2
    dpi = 1200;
end

if ispc
    downloads_dir = fullfile(getenv('USERPROFILE'),'Downloads');
else
    downloads_dir = fullfile(getenv('HOME'),'Downloads');
end
if ~exist(downloads_dir,'dir')
    mkdir(downloads_dir);
end

outpath = fullfile(downloads_dir, fname);

try
    exportgraphics(gcf, outpath, 'Resolution', dpi);
    fprintf('Saved %s to %s\n', fname, downloads_dir);
catch
    warning('exportgraphics failed. Falling back to print.');
    print(gcf, outpath, '-dpng', ['-r' num2str(dpi)]);
end
end

function [saved_states,saved_logprob] = DEMCz( ...
    logtarget, states, Z, n, n_thin, blockindex, ...
    current_logprob, temperature_schedule, gamma_schedule, jitter)

[m,d] = size(states);
mZ0 = size(Z,1);

if nargin < 10 || isempty(jitter)
    jitter = 1e-5 .* ones(1,d);
end

if nargin < 9 || isempty(gamma_schedule)
    gamma_schedule = zeros(numel(blockindex),1);
    for k = 1:numel(blockindex)
        gamma_schedule(k) = 2.38 / sqrt(2 * numel(blockindex{k}));
    end
end

if isvector(gamma_schedule)
    gamma_schedule = gamma_schedule(:);
end

if nargin < 8 || isempty(temperature_schedule)
    temperature_schedule = [];
end
if nargin < 7 || isempty(current_logprob)
    current_logprob = zeros(m,1);
    for j = 1:m
        current_logprob(j) = logtarget(states(j,:));
    end
end
if nargin < 6 || isempty(blockindex)
    blockindex = {1:d};
end
if nargin < 5 || isempty(n_thin)
    n_thin = 1;
end

temperature_schedule = temperature_schedule(:);
if numel(temperature_schedule) < n
    temperature_schedule = [temperature_schedule; ones(n-numel(temperature_schedule),1)];
end

gamidx = mod(0:(n-1), size(gamma_schedule,2)) + 1;
save_logprob = [];

for i = 1:n
    for j = 1:m
        for k = 1:numel(blockindex)
            idx = blockindex{k};
            mZ = size(Z,1) - m;

            if mZ < 2
                n1 = randi(size(Z,1));
                if size(Z,1) < 2
                    n2 = n1;
                else
                    n2 = randi(size(Z,1)-1);
                    if n2 >= n1
                        n2 = n2 + 1;
                    end
                end
            else
                n1 = randi(mZ);
                n2 = randi(mZ-1);
                if n2 >= n1
                    n2 = n2 + 1;
                end
            end

            diff = zeros(1,d);
            gamma_k = gamma_schedule(k, gamidx(i));
            diff(idx) = gamma_k .* (Z(n1,idx) - Z(n2,idx)) + ...
                        jitter(idx) .* randn(1,numel(idx));

            proposal = states(j,:) + diff;
            prop_logp = logtarget(proposal);

            if isfinite(prop_logp) && ...
               (temperature_schedule(i) * log(rand) < (prop_logp - current_logprob(j)))
                states(j,:) = proposal;
                current_logprob(j) = prop_logp;
            end
        end
    end

    if mod(i, n_thin) == 0
        Z = [Z; states];
        save_logprob = [save_logprob; current_logprob];
    end
end

num_discard = mZ0;
if size(Z,1) > num_discard
    saved_states = Z((num_discard+1):end,:);
else
    saved_states = [];
end

saved_logprob = save_logprob;
end