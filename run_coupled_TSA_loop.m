%% run_coupled_TSA_loop.m
% Coupled adsorber-desorber TSA loop
% First simple version: fixed sorbent flow and regeneration gas flow
% Goal: connect beta_lean -> adsorber -> beta_rich -> desorber -> beta_lean

clear; clc;

%% =========================
% FIXED INPUTS
% ==========================

% Number of stages
N_ads = 4;
N_des = 4;

% Sorbent / isotherm model
Isotherm = 'Toth_Lewatit_Zerobin';
%Isotherm = 'Toth_Lewatit_Veneman';
%Isotherm = 'Toth_Lewatit_Sutanto';
%Isotherm = 'Toth_13X_Zerobin';

% Operating temperatures
T_ads = 50 + 273.15;      % K
T_des = 120 + 273.15;     % K

% Pressure
P_bar = 1.01325;

% Gas feed to adsorber
n_g_feed = 1.0;      % mol/s total flue gas feed
yCO2_feed = 0.10;    % inlet CO2 mole fraction

% CO2 capture target
capture_target = 0.90;

% Molecular weights
MW_CO2 = 44.0095e-3;    % kg/mol
MW_H2O = 18.01528e-3;   % kg/mol

% Target captured CO2 flow
n_CO2_target = capture_target * n_g_feed * yCO2_feed;   % mol/s
m_CO2_target = n_CO2_target * MW_CO2;                   % kg/s
% Sorbent circulation rate
m_s = 0.04914;          % kg/s

% Stripping steam to desorber
% Steam is treated as a non-adsorbing gas because H2O adsorption
% is neglected in the current equilibrium model.

% Temporary validation value:
% selected to reproduce the previous n_g_regen = 0.20 mol/s case
steam_to_CO2_ratio = 0.9096662211;   % kg steam/kg CO2 captured

m_steam_in = steam_to_CO2_ratio * m_CO2_target;   % kg/s
n_steam_in = m_steam_in / MW_H2O;                 % mol/s

yCO2_steam_in = 0.0;   % pure stripping steam

%% =========================
% LOOP SETTINGS
% ==========================

% Supervisor advice: start with very low lean loading
beta_lean_guess = 1e-4;      % mol/kg

max_iter = 50;
tol_beta = 1e-6;

% Optional relaxation to avoid oscillation
relaxation = 0.5;

%% =========================
% ITERATIVE COUPLED LOOP
% ==========================

beta_lean_old = beta_lean_guess;

fprintf('\n===== COUPLED TSA LOOP ITERATION =====\n');

for iter = 1:max_iter

    %% -------------------------
    % 1. Run adsorber
    % --------------------------

    ads_results = solve_adsorber_Nstage(N_ads, n_g_feed, yCO2_feed, ...
                                        m_s, beta_lean_old, ...
                                        T_ads, P_bar, ...
                                        Isotherm);

    beta_rich = ads_results.beta_rich;

    %% -------------------------
    % 2. Run desorber
    % --------------------------
    
    des_results = solve_desorber_Nstage(N_des, n_steam_in, yCO2_steam_in, ...
                                    m_s, beta_rich, ...
                                    T_des, P_bar, ...
                                    Isotherm);

    beta_lean_calculated = des_results.beta_lean;

    %% -------------------------
    % 3. Check convergence
    % --------------------------

    difference = beta_lean_calculated - beta_lean_old;

    fprintf('Iter %2d: beta_lean_old = %.8f, beta_rich = %.8f, beta_lean_new = %.8f, diff = %.3e, capture = %.2f %%\n', ...
            iter, beta_lean_old, beta_rich, beta_lean_calculated, difference, ...
            100 * ads_results.capture_fraction);

    if abs(difference) < tol_beta
        fprintf('\nLoop converged after %d iterations.\n', iter);
        break;
    end

    %% -------------------------
    % 4. Update lean loading
    % --------------------------
    % Relaxed update:
    % beta_lean_old is moved partially toward beta_lean_calculated

    beta_lean_old = beta_lean_old + relaxation * difference;

end

if iter == max_iter
    fprintf('\nWarning: loop did not converge within max_iter.\n');
end

%% =========================
% FINAL RESULTS
% ==========================

fprintf('\n===== FINAL COUPLED TSA RESULTS =====\n');

fprintf('Isotherm                         = %s\n', Isotherm);
fprintf('Adsorber stages                  = %d\n', N_ads);
fprintf('Desorber stages                  = %d\n', N_des);
fprintf('Adsorber temperature             = %.2f K\n', T_ads);
fprintf('Desorber temperature             = %.2f K\n', T_des);
fprintf('Sorbent circulation rate         = %.6f kg/s\n', m_s);

fprintf('\nLoadings:\n');
fprintf('Lean loading to adsorber         = %.8f mol/kg\n', des_results.beta_lean);
fprintf('Rich loading from adsorber       = %.8f mol/kg\n', ads_results.beta_rich);
fprintf('Working capacity                 = %.8f mol/kg\n', ads_results.beta_rich - des_results.beta_lean);

fprintf('\nAdsorber performance:\n');
fprintf('Feed yCO2                        = %.6f\n', yCO2_feed);
fprintf('Outlet yCO2                      = %.6f\n', ads_results.y_gas_out);
fprintf('Capture fraction                 = %.2f %%\n', 100 * ads_results.capture_fraction);
fprintf('CO2 captured                     = %.8f mol/s\n', ads_results.n_CO2_captured);
fprintf('Adsorber cooling duty            = %.2f W\n', ads_results.Q_ads);

fprintf('\nDesorber performance:\n');
fprintf('Stripping steam flow             = %.8f mol/s\n', n_steam_in);
fprintf('Stripping steam mass flow        = %.8f kg/s\n', m_steam_in);
fprintf('Steam-to-CO2 stripping ratio     = %.8f kg steam/kg CO2\n', ...
    steam_to_CO2_ratio);
fprintf('Steam inlet yCO2                 = %.6f\n', yCO2_steam_in);
fprintf('CO2-rich gas yCO2 out            = %.6f\n', des_results.y_CO2_rich_gas_out);
fprintf('CO2 desorbed                     = %.8f mol/s\n', des_results.n_CO2_desorbed);
fprintf('Desorber heat duty               = %.2f W\n', des_results.Q_des);
%% =========================
% MASS BALANCE CHECK
% ==========================

CO2_by_solid_ads = m_s * ...
    (ads_results.beta_rich - ads_results.beta_lean);

CO2_by_solid_des = m_s * ...
    (des_results.beta_rich - des_results.beta_lean);

fprintf('\nMass balance check:\n');

fprintf('CO2 captured from gas in adsorber = %.8f mol/s\n', ...
    ads_results.n_CO2_captured);

fprintf('CO2 uptake by solid in adsorber   = %.8f mol/s\n', ...
    CO2_by_solid_ads);

fprintf('CO2 desorbed to gas in desorber   = %.8f mol/s\n', ...
    des_results.n_CO2_desorbed);

fprintf('CO2 released by solid in desorber = %.8f mol/s\n', ...
    CO2_by_solid_des);

fprintf('Adsorber gas-solid mismatch       = %.8e mol/s\n', ...
    ads_results.n_CO2_captured - CO2_by_solid_ads);

fprintf('Loop capture-desorption mismatch  = %.8e mol/s\n', ...
    ads_results.n_CO2_captured - des_results.n_CO2_desorbed);
%% =========================
% RECIRCULATION RATIOS
% ==========================

R_s_gas = m_s / n_g_feed;
R_s_CO2 = m_s / (n_g_feed * yCO2_feed);

fprintf('\nRecirculation ratios:\n');
fprintf('Sorbent-to-gas ratio             = %.8f kg sorbent / mol gas\n', R_s_gas);
fprintf('Sorbent-to-inlet-CO2 ratio       = %.8f kg sorbent / mol CO2,in\n', R_s_CO2);

%% =========================
% BASIC ENERGY INDICATOR
% ==========================

Q_basic_total = ads_results.Q_ads + des_results.Q_des;

if ads_results.m_CO2_captured > 0
    q_basic_specific = 1e-6 * Q_basic_total / ads_results.m_CO2_captured;
else
    q_basic_specific = NaN;
end

fprintf('\nBasic energy indicator:\n');
fprintf('Q_ads + Q_des                    = %.2f W\n', Q_basic_total);
fprintf('Specific basic energy demand     = %.4f MJ/kg CO2\n', q_basic_specific);

%% =========================
% PLOTS
% ==========================

stage_ads = 1:N_ads;
stage_des = 1:N_des;

figure;
plot(stage_ads, ads_results.y_profile, 'o-', 'LineWidth', 1.5);
xlabel('Adsorber stage number');
ylabel('Gas-phase y_{CO2}');
title('Coupled TSA: Adsorber Gas CO2 Profile');
grid on;

figure;
plot(stage_ads, ads_results.beta_profile, 's-', 'LineWidth', 1.5);
xlabel('Adsorber stage number');
ylabel('\beta_{CO2} [mol/kg]');
title('Coupled TSA: Adsorber Sorbent Loading Profile');
grid on;

figure;
plot(stage_des, des_results.y_profile, 'o-', 'LineWidth', 1.5);
xlabel('Desorber stage number');
ylabel('Gas-phase y_{CO2}');
title('Coupled TSA: Desorber Gas CO2 Profile');
grid on;

figure;
plot(stage_des, des_results.beta_profile, 's-', 'LineWidth', 1.5);
xlabel('Desorber stage number');
ylabel('\beta_{CO2} [mol/kg]');
title('Coupled TSA: Desorber Sorbent Loading Profile');
grid on;

figure;
bar(stage_ads, ads_results.Q_ads_stage/1000);
xlabel('Adsorber stage number');
ylabel('Cooling duty [kW]');
title('Coupled TSA: Stagewise Adsorber Cooling Duty');
grid on;

figure;
bar(stage_des, des_results.Q_des_stage/1000);
xlabel('Desorber stage number');
ylabel('Heating duty [kW]');
title('Coupled TSA: Stagewise Desorber Heating Duty');
grid on;