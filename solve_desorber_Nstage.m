function results = solve_desorber_Nstage(N, n_steam_in, yCO2_steam_in, ...
                                         m_s, beta_rich, ...
                                         T_des, P_bar, ...
                                         Isotherm)
% solve_desorber_Nstage
% General N-stage counter-current equilibrium desorber model.
%
% Stage indexing:
%   Solid flows: Stage 1 -> Stage 2 -> ... -> Stage N
%   Gas flows:   Stage N -> Stage N-1 -> ... -> Stage 1
%
% Therefore:
%   Rich sorbent enters Stage 1
%   Lean sorbent leaves Stage N
%   Stripping steam enters Stage N
%   CO2-rich gas leaves Stage 1

   %% =========================
% INITIAL GUESSES
% ==========================

% Gas CO2 is expected to be highest at Stage 1
% and lowest close to the regeneration gas inlet at Stage N.
y_high_guess = min(0.80, ...
                   max(0.05, beta_rich * m_s / ...
                   (n_steam_in + beta_rich * m_s)));

y_low_guess = max(yCO2_steam_in * 10, 1e-5);

y_guess = linspace(y_high_guess, y_low_guess, N);

% Calculate physically meaningful loading guesses
% from the selected equilibrium isotherm.
beta_guess = zeros(1, N);

for i = 1:N
    pCO2_guess_i = y_guess(i) * P_bar;

    beta_guess(i) = feval(Isotherm, ...
                          pCO2_guess_i, ...
                          T_des);
end

% Unknown vector:
% x(1:N)       = gas outlet yCO2 from each stage
% x(N+1:2N)    = solid outlet beta from each stage
x0 = [y_guess, beta_guess];
    %% =========================
    % SOLVE NONLINEAR SYSTEM
    % ==========================

    options = optimset('Display', 'iter', ...
                       'TolFun', 1e-10, ...
                       'TolX', 1e-10, ...
                       'MaxIter', 1000, ...
                       'MaxFunEvals', 10000);

    residual_fun = @(x) desorber_Nstage_residual(x, N, ...
                                                  n_steam_in, yCO2_steam_in, ...
                                                  m_s, beta_rich, ...
                                                  T_des, P_bar, ...
                                                  Isotherm);

    [x_sol, fval, exitflag] = fsolve(residual_fun, x0, options);
    imag_x_max = max(abs(imag(x_sol)));
    imag_f_max = max(abs(imag(fval)));

    fprintf('\nDesorber complex-number check:\n');
    fprintf('Maximum imaginary part of x_sol = %.8e\n', imag_x_max);
    fprintf('Maximum imaginary part of fval  = %.8e\n', imag_f_max);

    if imag_x_max > 1e-10 || imag_f_max > 1e-10
        error(['Desorber solver returned a significant complex solution. ', ...
            'max imag(x) = %.6e, max imag(fval) = %.6e'], ...
            imag_x_max, imag_f_max);
    end

    x_sol = real(x_sol);
    fval = real(fval);
fprintf('\n===== DESORBER FSOLVE DIAGNOSTIC =====\n');
fprintf('Exit flag                   = %d\n', exitflag);
fprintf('Residual norm               = %.8e\n', norm(fval));
fprintf('Maximum absolute residual   = %.8e\n', max(abs(fval)));

if exitflag <= 0
    error(['Desorber fsolve did not converge. ', ...
           'Exit flag = %d, residual norm = %.6e'], ...
           exitflag, norm(fval));
end

if max(abs(fval)) > 1e-7
    error(['Desorber solution does not satisfy the residual equations. ', ...
           'Maximum residual = %.6e'], ...
           max(abs(fval)));
end

    %% =========================
    % EXTRACT SOLUTION
    % ==========================

    y = x_sol(1:N);
    beta = x_sol(N+1:2*N);

    %% =========================
    % GAS FLOW PROFILE AND STAGEWISE CO2 DESORPTION
    % ==========================

    n_g_out = zeros(1, N);

    n_CO2_desorbed_stage = zeros(1, N);
    m_CO2_desorbed_stage = zeros(1, N);
    Q_des_stage = zeros(1, N);

    % Because gas flows from Stage N to Stage 1, calculate from N down to 1
    for i = N:-1:1

        % Gas inlet to stage i
        if i == N
            y_in = yCO2_steam_in;
            n_in = n_steam_in;
        else
            y_in = y(i+1);
            n_in = n_g_out(i+1);
        end

        % Gas outlet from stage i
        y_out = y(i);

        % Outlet gas flow from H2O component balance
        n_out = n_in * (1 - y_in) / (1 - y_out);
        n_g_out(i) = n_out;

        % CO2 entering and leaving gas phase in this stage
        n_CO2_in_i = n_in * y_in;
        n_CO2_out_i = n_out * y_out;

        % CO2 desorbed into gas in this stage
        n_CO2_desorbed_stage(i) = n_CO2_out_i - n_CO2_in_i;

        % CO2 desorbed in kg/s
        m_CO2_desorbed_stage(i) = n_CO2_desorbed_stage(i) * 44.0095e-3;

        % Stagewise desorption heat duty
        % Solid inlet loading to stage i
if i == 1
    beta_in_i = beta_rich;
else
    beta_in_i = beta(i-1);
end

% Solid outlet loading from stage i
beta_out_i = beta(i);

% Loading-dependent heat of desorption
% Same equilibrium heat magnitude, evaluated at desorber temperature.
heat_fun = @(beta_local) ...
    Toth_heat_adsorption(beta_local, T_des, Isotherm);

% Stagewise desorption heating duty
Q_des_stage(i) = m_s * integral( ...
    heat_fun, ...
    beta_out_i, ...
    beta_in_i);
    end

    %% =========================
    % OVERALL PERFORMANCE
    % ==========================

    y_CO2_rich_gas_out = y(1);
    n_g_rich_gas_out = n_g_out(1);

    beta_lean = beta(N);
    working_capacity = beta_rich - beta_lean;

    n_CO2_desorbed = sum(n_CO2_desorbed_stage);
    m_CO2_desorbed = n_CO2_desorbed * 44.0095e-3;

    Q_des = sum(Q_des_stage);

    if m_CO2_desorbed > 0
        q_des_specific = 1e-6 * Q_des / m_CO2_desorbed;   % MJ/kg CO2
    else
        q_des_specific = NaN;
    end

    %% =========================
    % STORE RESULTS
    % ==========================

    results.N = N;

    results.y_profile = y;
    results.beta_profile = beta;
    results.n_g_profile = n_g_out;

    results.y_CO2_rich_gas_out = y_CO2_rich_gas_out;
    results.n_g_rich_gas_out = n_g_rich_gas_out;

    results.beta_rich = beta_rich;
    results.beta_lean = beta_lean;
    results.working_capacity = working_capacity;

    results.n_CO2_desorbed = n_CO2_desorbed;
    results.m_CO2_desorbed = m_CO2_desorbed;

    results.n_CO2_desorbed_stage = n_CO2_desorbed_stage;
    results.m_CO2_desorbed_stage = m_CO2_desorbed_stage;

    results.Q_des_stage = Q_des_stage;
    results.Q_des = Q_des;
    results.q_des_specific = q_des_specific;
end
