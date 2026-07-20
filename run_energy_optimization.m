clear; clc;

%% =========================
% ENERGY OPTIMIZATION SCRIPT
% ==========================

%% Model settings

Isotherm = 'Toth_Lewatit_Veneman';

N_ads = 4;
N_des = 4;

T_ads = 323.15;      % K, 50 C
T_des = 393.15;      % K, 120 C
P_bar = 1.01325;     % bar

% Adsorber feed, dry basis
n_g_feed = 1.0;      % mol/s dry gas
yCO2_feed = 0.10;    % dry CO2 mole fraction

% Capture target
capture_target = 0.90;

%% Molecular weights

MW_CO2 = 44.01e-3;      % kg/mol
MW_H2O = 18.01528e-3;   % kg/mol

%% Target captured CO2 mass flow

n_CO2_feed = n_g_feed * yCO2_feed;
n_CO2_target = capture_target * n_CO2_feed;
m_CO2_target = n_CO2_target * MW_CO2;

%% Coupled-loop settings

beta_lean_guess = 1e-4;
relaxation = 0.3;
tol_beta = 1e-6;
max_iter = 100;

%% Sorbent circulation search bracket

m_s_bracket = [0.03, 0.18];   % kg/s

%% Steam-to-CO2 ratios to test

steam_to_CO2_ratio_list = [0.02 0.03 0.04 0.05 0.06 0.07 0.08 ...
                           0.10 0.12 0.15 0.20 0.25 0.30 0.40 0.50];
% kg steam / kg captured CO2

%% Storage

n_cases = length(steam_to_CO2_ratio_list);

results_table = table( ...
    zeros(n_cases,1), ...
    zeros(n_cases,1), ...
    zeros(n_cases,1), ...
    zeros(n_cases,1), ...
    zeros(n_cases,1), ...
    zeros(n_cases,1), ...
    zeros(n_cases,1), ...
    zeros(n_cases,1), ...
    zeros(n_cases,1), ...
    zeros(n_cases,1), ...
    'VariableNames', { ...
        'steam_to_CO2_ratio', ...
        'm_s_kg_s', ...
        'specific_solid_circulation', ...
        'capture_fraction', ...
        'beta_lean', ...
        'beta_rich', ...
        'working_capacity', ...
        'q_desorption', ...
        'q_sensible_net', ...
        'q_steam'});

results_table.q_total = zeros(n_cases,1);

%% =========================
% OPTIMIZATION LOOP
% ==========================

for k = 1:n_cases

    steam_to_CO2_ratio = steam_to_CO2_ratio_list(k);

    m_steam_in = steam_to_CO2_ratio * m_CO2_target;
    n_steam_in = m_steam_in / MW_H2O;

    fprintf('\n=============================================\n');
    fprintf('Steam-to-CO2 ratio = %.4f kg/kg\n', steam_to_CO2_ratio);
    fprintf('=============================================\n');

    try

        results = target_sorbent_circulation( ...
            N_ads, N_des, ...
            n_g_feed, yCO2_feed, ...
            n_steam_in, ...
            capture_target, ...
            m_s_bracket, ...
            T_ads, T_des, P_bar, ...
            Isotherm, ...
            beta_lean_guess, ...
            relaxation, tol_beta, max_iter);

        results_table.steam_to_CO2_ratio(k) = steam_to_CO2_ratio;
        results_table.m_s_kg_s(k) = results.m_s;
        results_table.specific_solid_circulation(k) = ...
            results.m_s / m_CO2_target;

        results_table.capture_fraction(k) = results.capture_fraction;
        results_table.beta_lean(k) = results.beta_lean;
        results_table.beta_rich(k) = results.beta_rich;
        results_table.working_capacity(k) = results.working_capacity;

        results_table.q_desorption(k) = results.q_desorption_specific;
        results_table.q_sensible_net(k) = results.q_sorbent_sensible_net_specific;
        results_table.q_steam(k) = results.q_steam_generation_specific;
        results_table.q_total(k) = results.q_total_heat_specific;

    catch ME

        warning('Case failed for steam ratio %.4f: %s', ...
            steam_to_CO2_ratio, ME.message);

        results_table.steam_to_CO2_ratio(k) = steam_to_CO2_ratio;
        results_table{k, 2:end} = NaN;

    end

end

%% =========================
% FIND MINIMUM ENERGY DEMAND
% ==========================

valid_rows = ~isnan(results_table.q_total);

valid_results = results_table(valid_rows, :);

[minimum_energy, idx_min] = min(valid_results.q_total);

best_case = valid_results(idx_min, :);

%% =========================
% PRINT RESULTS
% ==========================

fprintf('\n\n===== ENERGY OPTIMIZATION RESULTS =====\n');
disp(results_table);

fprintf('\n===== MINIMUM ENERGY DEMAND CASE =====\n');
disp(best_case);

fprintf('Minimum specific total heat demand = %.4f MJ/kg CO2\n', ...
    minimum_energy);

%% =========================
% PLOT 1: ENERGY VS STEAM RATIO
% ==========================

figure;
plot(valid_results.steam_to_CO2_ratio, ...
     valid_results.q_total, '-o', ...
     'LineWidth', 1.5);
hold on;

plot(best_case.steam_to_CO2_ratio, ...
     best_case.q_total, 'o', ...
     'MarkerSize', 10, ...
     'LineWidth', 2);

xlabel('Steam-to-CO2 ratio [kg steam/kg CO2]');
ylabel('Specific total heat demand [MJ/kg CO2]');
title('Energy demand vs steam-to-CO2 ratio for 90% CO2 capture');
grid on;

text(best_case.steam_to_CO2_ratio, ...
     best_case.q_total, ...
     sprintf('  Minimum = %.4f MJ/kg CO2', best_case.q_total));


%% =========================
% PLOT 2: ENERGY VS SOLID CIRCULATION
% ==========================

figure;
plot(valid_results.specific_solid_circulation, ...
     valid_results.q_total, '-o', ...
     'LineWidth', 1.5);
hold on;

plot(best_case.specific_solid_circulation, ...
     best_case.q_total, 'o', ...
     'MarkerSize', 10, ...
     'LineWidth', 2);

xlabel('Specific solid circulation [kg sorbent/kg CO2]');
ylabel('Specific total heat demand [MJ/kg CO2]');
title('Energy demand vs solid circulation for 90% CO2 capture');
grid on;

text(best_case.specific_solid_circulation, ...
     best_case.q_total, ...
     sprintf('  Minimum = %.4f MJ/kg CO2', best_case.q_total));
%% =========================
% PLOT 3: 3D OPTIMIZATION CURVE
% ==========================

figure;

plot3(valid_results.steam_to_CO2_ratio, ...
    valid_results.specific_solid_circulation, ...
    valid_results.q_total, '-o', ...
    'LineWidth', 1.5);

hold on;

plot3(best_case.steam_to_CO2_ratio, ...
    best_case.specific_solid_circulation, ...
    best_case.q_total, 'o', ...
    'MarkerSize', 10, ...
    'LineWidth', 2);

xlabel('Steam-to-CO2 ratio [kg steam/kg CO2]');
ylabel('Specific solid circulation [kg sorbent/kg CO2]');
zlabel('Specific total heat demand [MJ/kg CO2]');

title('Energy demand as function of steam ratio and solid circulation');
grid on;
view(45, 25);

text(best_case.steam_to_CO2_ratio, ...
    best_case.specific_solid_circulation, ...
    best_case.q_total, ...
    sprintf('  Minimum = %.4f MJ/kg CO2', best_case.q_total));
