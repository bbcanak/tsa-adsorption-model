function results = solve_coupled_TSA_fixed_ms( ...
    N_ads, N_des, ...
    n_g_feed, yCO2_feed, ...
    m_s, ...
    n_steam_in, ...
    T_ads, T_des, P_bar, ...
    Isotherm, ...
    beta_lean_guess, ...
    relaxation, tol_beta, max_iter)
% solve_coupled_TSA_fixed_ms
% Solves the closed adsorber-desorber TSA loop for a fixed sorbent
% circulation rate and fixed stripping steam flow.
%
% Coupled solution sequence:
%   beta_lean guess
%       -> adsorber
%       -> beta_rich
%       -> desorber
%       -> beta_lean calculated
%       -> relaxed beta_lean update
%       -> repeat until convergence
%
% The adsorber and desorber stage models are not modified here.

    %% =========================
    % FIXED STEAM INLET COMPOSITION
    % ==========================

    yCO2_steam_in = 0.0;

    %% =========================
    % INITIALIZE CLOSED LOOP
    % ==========================

    beta_lean_old = beta_lean_guess;
    converged = false;

    %% =========================
    % ITERATIVE COUPLED LOOP
    % ==========================

    for iter = 1:max_iter

        % -------------------------
        % 1. Solve adsorber
        % -------------------------

        ads_results = solve_adsorber_Nstage( ...
            N_ads, ...
            n_g_feed, yCO2_feed, ...
            m_s, beta_lean_old, ...
            T_ads, P_bar, ...
            Isotherm);

        beta_rich = ads_results.beta_rich;

        % -------------------------
        % 2. Solve desorber
        % -------------------------

        des_results = solve_desorber_Nstage( ...
            N_des, ...
            n_steam_in, yCO2_steam_in, ...
            m_s, beta_rich, ...
            T_des, P_bar, ...
            Isotherm);

        beta_lean_calculated = des_results.beta_lean;

        % -------------------------
        % 3. Check convergence
        % -------------------------

        difference = beta_lean_calculated - beta_lean_old;

        if abs(difference) < tol_beta
            converged = true;
            break;
        end

        % -------------------------
        % 4. Relaxed lean-loading update
        % -------------------------

        beta_lean_old = beta_lean_old + ...
            relaxation * difference;

    end

    %% =========================
    % CONVERGENCE CHECK
    % ==========================

    if ~converged
        error(['Coupled TSA loop did not converge within %d iterations. ', ...
               'Final beta difference = %.6e mol/kg'], ...
               max_iter, difference);
    end

    %% =========================
    % LOOP MASS BALANCE
    % ==========================

    loop_CO2_mismatch = ...
        ads_results.n_CO2_captured ...
        - des_results.n_CO2_desorbed;

    %% =========================
    % STORE RESULTS
    % ==========================

    results.ads = ads_results;
    results.des = des_results;

    results.beta_lean = des_results.beta_lean;
    results.beta_rich = ads_results.beta_rich;

    results.working_capacity = ...
        results.beta_rich - results.beta_lean;

    results.capture_fraction = ...
        ads_results.capture_fraction;

    results.n_CO2_captured = ...
        ads_results.n_CO2_captured;

    results.m_CO2_captured = ...
        ads_results.m_CO2_captured;

    results.n_CO2_desorbed = ...
        des_results.n_CO2_desorbed;

    results.loop_CO2_mismatch = ...
        loop_CO2_mismatch;

    results.m_s = m_s;
    results.n_steam_in = n_steam_in;

    results.iterations = iter;
    results.beta_difference = difference;
    results.converged = converged;

end