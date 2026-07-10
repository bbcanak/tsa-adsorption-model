function results = solve_adsorber_Nstage(N, n_g_feed, yCO2_feed, ...
                                         m_s, beta_lean, ...
                                         T_ads, P_bar, ...
                                         Isotherm)
% solve_adsorber_Nstage
% General N-stage counter-current equilibrium adsorber model.
%
% Gas direction:
%   Stage 1 -> Stage 2 -> ... -> Stage N
%
% Solid direction:
%   Stage N -> Stage N-1 -> ... -> Stage 1
%
% Inputs:
%   N             number of equilibrium stages
%   n_g_feed      inlet gas flow, mol/s
%   yCO2_feed     inlet gas CO2 mole fraction
%   m_s           sorbent mass flow, kg/s
%   beta_lean     CO2 loading of lean sorbent entering stage N, mol/kg
%   T_ads         adsorber temperature, K
%   P_bar         pressure, bar
%   Isotherm      isotherm function name as string
%
% Outputs:
%   results       structure containing profiles and performance indicators

   %% =========================
% INITIAL GUESSES
% ==========================

% Initial gas-phase CO2 profile
% Gas flows from Stage 1 to Stage N and loses CO2.
y_guess = linspace(0.8 * yCO2_feed, ...
                   max(0.1 * yCO2_feed, 1e-5), N);

% Initial sorbent-loading guess based directly on the selected isotherm.
% This gives fsolve a physically meaningful equilibrium starting point.
beta_guess = zeros(1, N);

for i = 1:N
    pCO2_guess_i = y_guess(i) * P_bar;

    beta_guess(i) = feval(Isotherm, ...
                          pCO2_guess_i, ...
                          T_ads);
end

% Unknown vector:
% x(1:N)       = gas outlet yCO2 from each stage
% x(N+1:2*N)   = solid outlet beta from each stage
x0 = [y_guess, beta_guess];

    %% =========================
    % SOLVE NONLINEAR SYSTEM
    % ==========================

    options = optimset('Display', 'iter', ...
                       'TolFun', 1e-10, ...
                       'TolX', 1e-10, ...
                       'MaxIter', 1000, ...
                       'MaxFunEvals', 10000);

    residual_fun = @(x) adsorber_Nstage_residual(x, N, ...
                                                  n_g_feed, yCO2_feed, ...
                                                  m_s, beta_lean, ...
                                                  T_ads, P_bar, ...
                                                  Isotherm);

    [x_sol, fval, exitflag] = fsolve(residual_fun, x0, options);
    imag_x_max = max(abs(imag(x_sol)));
    imag_f_max = max(abs(imag(fval)));

    fprintf('\nAdsorber complex-number check:\n');
    fprintf('Maximum imaginary part of x_sol = %.8e\n', imag_x_max);
    fprintf('Maximum imaginary part of fval  = %.8e\n', imag_f_max);

    if imag_x_max > 1e-10 || imag_f_max > 1e-10
        error(['Adsorber solver returned a significant complex solution. ', ...
            'max imag(x) = %.6e, max imag(fval) = %.6e'], ...
            imag_x_max, imag_f_max);
    end

    x_sol = real(x_sol);
    fval = real(fval);
fprintf('\n===== FSOLVE DIAGNOSTIC =====\n');
fprintf('Exit flag                   = %d\n', exitflag);
fprintf('Residual norm               = %.8e\n', norm(fval));
fprintf('Maximum absolute residual   = %.8e\n', max(abs(fval)));
if exitflag <= 0
    error(['Adsorber fsolve did not converge. ', ...
        'Exit flag = %d, residual norm = %.6e'], ...
        exitflag, norm(fval));
end

if max(abs(fval)) > 1e-7
    error(['Adsorber solution does not satisfy the residual equations. ', ...
        'Maximum residual = %.6e'], ...
        max(abs(fval)));
end

    %% =========================
    % EXTRACT SOLUTION
    % ==========================

    y = x_sol(1:N);
    beta = x_sol(N+1:2*N);

%% =========================
% GAS FLOW PROFILE AND STAGEWISE CO2 CAPTURE
% ==========================

n_g = zeros(1, N);

n_CO2_captured_stage = zeros(1, N);
m_CO2_captured_stage = zeros(1, N);
Q_ads_stage = zeros(1, N);

for i = 1:N

    % Gas inlet to stage i
    if i == 1
        y_in = yCO2_feed;
        n_in = n_g_feed;
    else
        y_in = y(i-1);
        n_in = n_g(i-1);
    end

    % Gas outlet from stage i
    y_out = y(i);

    % Outlet gas flow from inert balance
    n_out = n_in * (1 - y_in) / (1 - y_out);
    n_g(i) = n_out;

    % CO2 entering and leaving this stage
    n_CO2_in_i = n_in * y_in;
    n_CO2_out_i = n_out * y_out;

    % CO2 captured in this stage
    n_CO2_captured_stage(i) = n_CO2_in_i - n_CO2_out_i;

    % CO2 captured in kg/s
    m_CO2_captured_stage(i) = ...
        n_CO2_captured_stage(i) * 44.0095e-3;

    % Solid inlet loading to stage i
    if i == N
        beta_in_i = beta_lean;
    else
        beta_in_i = beta(i+1);
    end

    % Solid outlet loading from stage i
    beta_out_i = beta(i);

    % Loading-dependent heat of adsorption
    heat_fun = @(beta_local) ...
        Toth_heat_adsorption(beta_local, T_ads, Isotherm);

    % Stagewise adsorption cooling duty
    Q_ads_stage(i) = m_s * integral( ...
        heat_fun, ...
        beta_in_i, ...
        beta_out_i);

end
    %% =========================
    % OVERALL PERFORMANCE
    % ==========================

    y_gas_out = y(N);
    n_g_out = n_g(N);

    n_CO2_in = n_g_feed * yCO2_feed;
    n_CO2_out = n_g_out * y_gas_out;

    n_CO2_captured = n_CO2_in - n_CO2_out;
    m_CO2_captured = n_CO2_captured * 44.0095e-3;   % kg/s

    capture_fraction = n_CO2_captured / n_CO2_in;

    beta_rich = beta(1);   % solid leaves adsorber from stage 1 as rich sorbent

    % Total adsorber cooling duty from stagewise sum
    Q_ads = sum(Q_ads_stage);   % W

    % Specific cooling duty
    q_ads_specific = 1e-6 * Q_ads / m_CO2_captured;   % MJ/kg CO2
    %% =========================
    % STORE RESULTS
    % ==========================

    results.N = N;
    results.y_profile = y;
    results.beta_profile = beta;
    results.n_g_profile = n_g;

    results.y_gas_out = y_gas_out;
    results.n_g_out = n_g_out;

    results.beta_lean = beta_lean;
    results.beta_rich = beta_rich;
    results.working_capacity = beta_rich - beta_lean;
    
    results.n_CO2_captured = n_CO2_captured;
    results.m_CO2_captured = m_CO2_captured;
    results.capture_fraction = capture_fraction;

    results.Q_ads = Q_ads;
    results.q_ads_specific = q_ads_specific;
    results.n_CO2_captured_stage = n_CO2_captured_stage;
    results.m_CO2_captured_stage = m_CO2_captured_stage;
    results.Q_ads_stage = Q_ads_stage;
    results.q_ads_stage_specific = 1e-6 * Q_ads_stage ./ m_CO2_captured_stage;
end