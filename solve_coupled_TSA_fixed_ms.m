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
% This function is used by the optimization routine.
% It returns capture performance and TSA heat-demand terms.

    %% =========================
    % CONSTANTS
    % ==========================

    MW_H2O = 18.01528e-3;   % kg/mol

    %% =========================
    % ENERGY MODEL SETTINGS
    % ==========================

    switch Isotherm

        case {'Toth_Lewatit_Zerobin', ...
              'Toth_Lewatit_Veneman', ...
              'Toth_Lewatit_Sutanto'}

            cp_sorbent = 1.58e3;      % J/kg/K

        case 'Toth_13X_Zerobin'

            cp_sorbent = 0.93e3;      % J/kg/K

        otherwise

            error('Unknown isotherm for sorbent heat capacity: %s', Isotherm);

    end

    heat_recovery_fraction = 0.50;

    T_feed_water = 100 + 273.15;       % K
    h_vap_H2O_100C = 2256.4e3;         % J/kg
    cp_H2O_vapor_ideal = 1.86e3;       % J/kg/K

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
    % TSA HEAT DEMAND
    % ==========================

    m_steam_in = n_steam_in * MW_H2O;

    Q_desorption = des_results.Q_des;

    Q_sorbent_sensible_gross = ...
        m_s * cp_sorbent * (T_des - T_ads);

    Q_sorbent_sensible_net = ...
        (1 - heat_recovery_fraction) * Q_sorbent_sensible_gross;

    Q_steam_generation = ...
        m_steam_in * ...
        (h_vap_H2O_100C + cp_H2O_vapor_ideal * (T_des - T_feed_water));

    Q_total_heat = ...
        Q_desorption ...
        + Q_sorbent_sensible_net ...
        + Q_steam_generation;

    if ads_results.m_CO2_captured > 0

        q_desorption_specific = ...
            1e-6 * Q_desorption / ads_results.m_CO2_captured;

        q_sorbent_sensible_gross_specific = ...
            1e-6 * Q_sorbent_sensible_gross / ads_results.m_CO2_captured;

        q_sorbent_sensible_net_specific = ...
            1e-6 * Q_sorbent_sensible_net / ads_results.m_CO2_captured;

        q_steam_generation_specific = ...
            1e-6 * Q_steam_generation / ads_results.m_CO2_captured;

        q_total_heat_specific = ...
            1e-6 * Q_total_heat / ads_results.m_CO2_captured;

    else

        q_desorption_specific = NaN;
        q_sorbent_sensible_gross_specific = NaN;
        q_sorbent_sensible_net_specific = NaN;
        q_steam_generation_specific = NaN;
        q_total_heat_specific = NaN;

    end

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
    results.m_steam_in = m_steam_in;

    results.cp_sorbent = cp_sorbent;
    results.heat_recovery_fraction = heat_recovery_fraction;

    results.Q_desorption = Q_desorption;
    results.Q_sorbent_sensible_gross = Q_sorbent_sensible_gross;
    results.Q_sorbent_sensible_net = Q_sorbent_sensible_net;
    results.Q_steam_generation = Q_steam_generation;
    results.Q_total_heat = Q_total_heat;

    results.q_desorption_specific = q_desorption_specific;
    results.q_sorbent_sensible_gross_specific = q_sorbent_sensible_gross_specific;
    results.q_sorbent_sensible_net_specific = q_sorbent_sensible_net_specific;
    results.q_steam_generation_specific = q_steam_generation_specific;
    results.q_total_heat_specific = q_total_heat_specific;

    results.iterations = iter;
    results.beta_difference = difference;
    results.converged = converged;
    % Store final stagewise adsorber/desorber results
    results.ads_results = ads_results;
    results.des_results = des_results;

end