function [Si_matrix, Ti_matrix] = sobolshi()


    % -------------------- Output folder: DOWNLOADS --------------------
    stamp = datestr(now,'yyyymmdd_HHMMSS');
    downloads_dir = get_downloads_dir();

    % Put everything in its own run folder to avoid overwrite/Excel locks
    outdir = fullfile(downloads_dir, ['SobolOutputs_' stamp]);
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end

    % -------------------- Time points (hours) --------------------
    time_points = [0, 4, 8, 12, 16, 24, 32, 40, 48, 56, 64, 72];
    time_points = unique(time_points(:))';  % ensure unique & row

    % -------------------- Parameter names (LaTeX) --------------------
    param_names = {'$k_1$', '$k_{-1}$', '$k_3$', '$k_{-3}$', '$k_{\mathrm{conv}}$', ...
                   '$k_{\mathrm{deg},C}$', '$k_{\mathrm{deg},A}$', '$D_{e0}$'};

    % -------------------- Bounds --------------------
    lb = [ ...
        5.0e-2, ...   % k1
        1.0e-1, ...   % k_{-1}
        1.0e-2, ...   % k3
        1.0e-3, ...   % k_{-3}
        5.0e-2, ...   % kconv
        1.0e-5, ...   % kdegC
        1.0e-4, ...   % kdegA
        5.0e-1  ...   % De0
    ];

    ub = [ ...
        2.0,    ...   % k1
        3.0,    ...   % k_{-1}
        1.0,    ...   % k3
        2.0e-1, ...   % k_{-3}
        2.0,    ...   % kconv
        1.0e-2, ...   % kdegC
        5.0e-2, ...   % kdegA
        2.0e1   ...   % De0
    ];

    lb = lb(:)'; ub = ub(:)';  % row vectors
    p  = numel(lb);
    assert(numel(ub) == p, 'lb and ub must have same length.');
    assert(all(ub > lb), 'All upper bounds must be > lower bounds.');

    % -------------------- Sobol settings --------------------
    sample_points = 256;    % M
    use_log1p = false;      % set true if P spans orders of magnitude

    settings.model         = @(par) simulate_P_PDE_mean(par, time_points, use_log1p);
    settings.lower_bounds  = lb;
    settings.upper_bounds  = ub;
    settings.sample_points = sample_points;

    fprintf('\n=== Sobol (covariance Si, subtract-from-1 Ti) for PDE Product P(t) ===\n');
    fprintf('Output folder (Downloads): %s\n', outdir);
    fprintf('Convention: C_j = B_A^(j) (start from B, replace column j with A)\n');
    fprintf('Si(t,j) = (E[f(A)f(Cj)] - f0^2)/Var(f)\n');
    fprintf('Ti(t,j) = 1 - (E[f(B)f(Cj)] - f0^2)/Var(f)\n');
    fprintf('Times (hours): %s\n', mat2str(time_points));
    fprintf('M=%d, p=%d => total PDE solves ≈ (p+2)*M = %d\n', ...
            sample_points, p, (p+2)*sample_points);

    % -------------------- Run Sobol --------------------
    out = saltelli_sobol_time_series_sub1(settings);
    Si_matrix = out.Si;   % Nt x p
    Ti_matrix = out.Ti;   % Nt x p

    % -------------------- Unique markers per parameter --------------------
    markers = {'o','s','d','^','v','>','<','p'};   % 8 unique markers
    marker_idx = 1:numel(time_points);            % set to 1:2:end for fewer markers

    % -------------------- Plot First-order indices --------------------
    fig1 = figure('Color','w','Position',[100 100 900 560]); hold on;
    for j = 1:p
        mk = markers{j};
        plot(time_points, Si_matrix(:,j), ['-' mk], ...
             'LineWidth', 2.5, 'MarkerSize', 6, 'MarkerIndices', marker_idx, ...
             'DisplayName', param_names{j});
    end
    xlabel('Time (hours)', 'FontSize', 16);
    ylabel('First-order Sobol Index', 'Interpreter','latex', 'FontSize', 16);
    ylim([0 1]); grid on; set(gca, 'FontSize', 14);
    lgd = legend(param_names, 'Interpreter','latex', 'FontSize', 11, 'NumColumns', 4);
    lgd.Location = 'northoutside'; lgd.Box = 'off';
    exportgraphics(fig1, fullfile(outdir, 'sobol_first_order_P.png'), 'Resolution', 600);

    % -------------------- Plot Total-order indices --------------------
    fig2 = figure('Color','w','Position',[100 100 900 560]); hold on;
    for j = 1:p
        mk = markers{j};
        plot(time_points, Ti_matrix(:,j), ['-' mk], ...
             'LineWidth', 2.5, 'MarkerSize', 6, 'MarkerIndices', marker_idx, ...
             'DisplayName', param_names{j});
    end
    xlabel('Time (hours)', 'FontSize', 16);
    ylabel('Total-order Sobol Index', 'Interpreter','latex', 'FontSize', 16);
    ylim([0 1]); grid on; set(gca, 'FontSize', 14);
    lgd = legend(param_names, 'Interpreter','latex', 'FontSize', 11, 'NumColumns', 4);
    lgd.Location = 'northoutside'; lgd.Box = 'off';
    exportgraphics(fig2, fullfile(outdir, 'sobol_total_order_P.png'), 'Resolution', 600);

    % -------------------- Save CSVs --------------------
    writematrix(time_points(:), fullfile(outdir, 'sobol_time_points_hours.csv'));
    writematrix(Si_matrix,      fullfile(outdir, 'sobol_Si_P.csv'));
    writematrix(Ti_matrix,      fullfile(outdir, 'sobol_Ti_P.csv'));

    % -------------------- Save Excel (.xlsx) --------------------
    excel_file = fullfile(outdir, 'sobol_indices_P.xlsx');

    varnames = {'Time_hours','k1','k_minus1','k3','k_minus3','kconv','kdegC','kdegA','De0'};
    TSi = array2table([time_points(:), Si_matrix], 'VariableNames', varnames);
    TTi = array2table([time_points(:), Ti_matrix], 'VariableNames', varnames);

    writetable(TSi, excel_file, 'Sheet', 'Si');
    writetable(TTi, excel_file, 'Sheet', 'Ti');

    param_rows = {'k1','k_minus1','k3','k_minus3','kconv','kdegC','kdegA','De0'}';
    Tbnds = table(lb(:), ub(:), 'VariableNames', {'lb','ub'}, 'RowNames', param_rows);
    writetable(Tbnds, excel_file, 'Sheet', 'bounds', 'WriteRowNames', true);

    Ttime = table(time_points(:), 'VariableNames', {'Time_hours'});
    writetable(Ttime, excel_file, 'Sheet', 'time_points');

    fprintf('✅ Saved everything to: %s\n', outdir);
end

% =========================================================================
% PDE simulation -> vector output at requested time points
% Output: spatial-mean Product P(x,t) (optionally log1p)
% =========================================================================
function P_vec = simulate_P_PDE_mean(pars, time_points, use_log1p)

    % --- MOL grid ---
    N_GRID   = 25;
    L_DOMAIN = 0.25;
    dx       = L_DOMAIN / (N_GRID - 1);

    % --- Initial/boundary constants ---
    C0_POLY = 0.302;
    A0_POLY = 0.698;
    E_INIT  = 0.0;
    P_INIT  = 0.0;

    E_BULK  = 1;
    ALPHA_E = 1.0;

    time_points = unique(time_points(:))';

    % States: [E, C, A, EC, EA, P] each length N_GRID
    E0  = E_INIT  * ones(1, N_GRID);
    C0  = C0_POLY * ones(1, N_GRID);
    A0  = A0_POLY * ones(1, N_GRID);
    EC0 = zeros(1, N_GRID);
    EA0 = zeros(1, N_GRID);
    P0  = P_INIT  * ones(1, N_GRID);

    E0(1) = E_BULK;

    y0   = [E0, C0, A0, EC0, EA0, P0]';
    opts = odeset('RelTol',1e-6,'AbsTol',1e-9,'MaxStep',0.5);

    rhs  = @(t,y) rhs_pde_mol(t, y, pars, N_GRID, dx, C0_POLY, A0_POLY, E_BULK, ALPHA_E);

    try
        [~, Y] = ode15s(rhs, time_points, y0, opts);
    catch
        P_vec = nan(1, numel(time_points));
        return;
    end

    idxP = (1:N_GRID) + 5*N_GRID;   % P is 6th block
    Pmat = Y(:, idxP);              % Nt x N_GRID

    P_vec = mean(Pmat, 2).';        % 1 x Nt

    if use_log1p
        P_vec = log1p(max(P_vec, 0));
    end
end

% =========================================================================
% PDE RHS (method-of-lines)
% =========================================================================
function dydt = rhs_pde_mol(~, y, pars, N_GRID, dx, C0_POLY, A0_POLY, E_BULK, ALPHA_E)

    k1       = pars(1);
    k_minus1 = pars(2);
    k3       = pars(3);
    k_minus3 = pars(4);
    kconv    = pars(5);
    kdegC    = pars(6);
    kdegA    = pars(7);
    De0      = pars(8);

    y  = y(:);

    E  = y(          1 :   N_GRID);
    C  = y(  N_GRID+1 : 2*N_GRID);
    A  = y(2*N_GRID+1 : 3*N_GRID);
    EC = y(3*N_GRID+1 : 4*N_GRID);
    EA = y(4*N_GRID+1 : 5*N_GRID);

    % Dirichlet boundary for enzyme at x=0
    E(1) = E_BULK;

    % Porosity-like factor (include bound states)
    denom0 = C0_POLY + A0_POLY;
    phi    = 1.0 - (C + A ) ./ max(denom0, 1e-12);
    phi    = min(max(phi, -0.5), 2.0);

    De    = max(De0 .* (1.0 + ALPHA_E .* phi), 1e-15);
    diffE = diffusion_dirichlet_neumann_varD(E, De, E_BULK, N_GRID, dx);

    dE  =  diffE ...
           - k1.*E.*C      + k_minus1.*EC ...
           - k3.*E.*A      + k_minus3.*EA ...
           + kconv.*EC ...
           + kdegC.*EC     + kdegA.*EA;

    dC  = -k1.*E.*C + k_minus1.*EC;
    dA  = -k3.*E.*A + k_minus3.*EA + kconv.*EC;
    dEC =  k1.*E.*C - k_minus1.*EC - kconv.*EC - kdegC.*EC;
    dEA =  k3.*E.*A - k_minus3.*EA - kdegA.*EA;

    dP  =  kdegC.*EC + kdegA.*EA;

    % Keep boundary fixed
    dE(1) = 0.0;

    dydt = [dE; dC; dA; dEC; dEA; dP];
end

% =========================================================================
% Variable-coefficient diffusion operator
% Dirichlet at left: u(1)=u_left
% Neumann at right: du/dx = 0
% =========================================================================
function div = diffusion_dirichlet_neumann_varD(u, D, u_left, N_GRID, dx)

    u    = u(:);
    D    = D(:);
    u(1) = u_left;

    flux           = zeros(N_GRID+1, 1);
    D_iface        = 0.5*(D(1:end-1) + D(2:end));
    flux(2:N_GRID) = D_iface .* (u(2:end) - u(1:end-1)) / dx;

    flux(N_GRID+1) = 0.0;  % Neumann at right

    div             = zeros(N_GRID, 1);
    div(2:N_GRID-1) = (flux(3:N_GRID)     - flux(2:N_GRID-1)) / dx;
    div(N_GRID)     = 2.0*(flux(N_GRID+1) - flux(N_GRID))     / dx;
    div(1)          = 0.0;
end

% =========================================================================
% Sobol indices (vector output) with total-order in "1 - ..." form
% Convention: C_j = B_A^(j)
% Returns out.Si (Nt x p), out.Ti (Nt x p)
% =========================================================================
function out = saltelli_sobol_time_series_sub1(settings)

    model = settings.model;
    lb    = settings.lower_bounds(:)';   % 1 x p
    ub    = settings.upper_bounds(:)';   % 1 x p
    M     = settings.sample_points;
    p     = numel(lb);

    rng(123);
    sob = sobolset(p, 'Skip', 1e3, 'Leap', 1e2);
    sob = scramble(sob, 'MatousekAffineOwen');
    AB  = net(sob, 2*M);

    A_mat = lb + AB(1:M,      :) .* (ub - lb);
    B_mat = lb + AB(M+1:2*M,  :) .* (ub - lb);

    % C_j = B with column j from A
    C_arr = zeros(M, p, p);
    for j = 1:p
        C_arr(:,:,j) = B_mat;
        C_arr(:,j,j) = A_mat(:,j);
    end

    y_test = model(A_mat(1,:));
    Nt = numel(y_test);

    yA = nan(M, Nt);
    yB = nan(M, Nt);
    yC = nan(M, p, Nt);

    for i = 1:M
        try, yA(i,:) = model(A_mat(i,:)); catch, end
        try, yB(i,:) = model(B_mat(i,:)); catch, end
    end

    for j = 1:p
        for i = 1:M
            try
                tmp = model(C_arr(i,:,j));
                yC(i,j,:) = reshape(tmp, 1, 1, Nt);
            catch
            end
        end
    end

    valid = all(isfinite(yA), 2) & all(isfinite(yB), 2);
    for j = 1:p
        valid = valid & all(isfinite(squeeze(yC(:,j,:))), 2);
    end

    n_valid = sum(valid);
    fprintf('Sobol valid samples: %d / %d\n', n_valid, M);

    if n_valid < max(50, floor(0.2*M))
        error('Too many failed PDE solves (%d valid of %d). Tighten bounds or reduce M.', n_valid, M);
    end

    yA = yA(valid, :);
    yB = yB(valid, :);
    yC = yC(valid, :, :);

    f0    = mean([yA; yB], 1);       % 1 x Nt
    var_y = var([yA; yB], 1, 1);     % 1 x Nt
    var_y = max(var_y, eps);

    Si = zeros(Nt, p);
    Ti = zeros(Nt, p);

    for j = 1:p
        yCj = squeeze(yC(:,j,:));    % N x Nt
        Si(:,j) = ( (mean(yA .* yCj, 1) - f0.^2) ./ var_y ).';
        Ti(:,j) = ( 1 - (mean(yB .* yCj, 1) - f0.^2) ./ var_y ).';
    end

    Si = max(min(Si, 1), 0);
    Ti = max(min(Ti, 1), 0);

    out.Si = Si;
    out.Ti = Ti;
end

% =========================================================================
% Downloads folder helper
% =========================================================================
function ddir = get_downloads_dir()
    if ispc
        ddir = fullfile(getenv('USERPROFILE'), 'Downloads');
    else
        ddir = fullfile(getenv('HOME'), 'Downloads');
    end
    if ~exist(ddir,'dir')
        mkdir(ddir);
    end
end