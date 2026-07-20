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
%Isotherm = 'Toth_Lewatit_Zerobin';
Isotherm = 'Toth_Lewatit_Veneman';
%Isotherm = 'Toth_Lewatit_Sutanto';
%Isotherm = 'Toth_13X_Zerobin';

% Operating temperatures
T_ads = 50 + 273.15;      % K
T_des = 120 + 273.15;     % K

% Pressure
P_bar = 1.01325;

% Gas feed to adsorber
% Current adsorber CO2/N2 model is kept on dry-gas basis.
n_g_feed = 1.0;      % mol/s dry flue gas feed
yCO2_feed = 0.10;    % dry inlet CO2 mole fraction

% Water vapor in adsorber feed for relative-humidity calculation
% H2O is treated as non-adsorbing.
yH2O_ads_feed_wet = 0.04;   % mol fraction H2O in wet adsorber feed
% CO2 capture target
capture_target = 0.90;

% Molecular weights
MW_CO2 = 44.0095e-3;    % kg/mol
MW_H2O = 18.01528e-3;   % kg/mol

% Target captured CO2 flow
n_CO2_target = capture_target * n_g_feed * yCO2_feed;   % mol/s
m_CO2_target = n_CO2_target * MW_CO2;                   % kg/s
%% =========================
% SORBENT CIRCULATION
% ==========================

% Zerobin Case 3 reported value:
% specific solids circulation rate =
% 19 kg unloaded sorbent / kg captured CO2

specific_solid_circulation = 19.0;     % kg sorbent / kg CO2 captured

m_s = specific_solid_circulation * m_CO2_target;   % kg/s

% Stripping steam to desorber
% Steam is treated as a non-adsorbing gas because H2O adsorption
% is neglected in the current equilibrium model.
% Zerobin Case 3 reported minimal stripping ratio:
% 0.08 kg steam per kg captured CO2.
% Here it is based on the target captured CO2 amount.

steam_to_CO2_ratio = 0.08;   % kg steam/kg CO2 captured

m_steam_in = steam_to_CO2_ratio * m_CO2_target;   % kg/s
n_steam_in = m_steam_in / MW_H2O;                 % mol/s

yCO2_steam_in = 0.0;   % pure stripping steam
%% =========================
% ENERGY MODEL SETTINGS
% ==========================

% Sorbent heat capacity
% Zerobin dissertation values:
% Lewatit VP OC 1065:
% RICH 1.68921 kJ/kg/K
% LEAN 1.9173 kJ/kg/K
% Zeolite 13X:        0.93 kJ/kg/K

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

% Lean-rich heat exchanger effectiveness
heat_recovery_fraction = 0.50;     % 50% heat recovery

% Feed-water and steam-generation assumptions
% First slim model following Tobias:
% liquid water at 100 C -> evaporation at 100 C -> ideal-gas H2O vapor at T_des
T_feed_water = 100 + 273.15;       % K

h_vap_H2O_100C = 2256.4e3;         % J/kg, latent heat at 100 C

cp_H2O_vapor_ideal = 1.86e3;       % J/kg/K, ideal-gas H2O vapor approx. near 100-120 C

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
% RELATIVE HUMIDITY PROFILES
% ==========================

T_ads_C = T_ads - 273.15;
T_des_C = T_des - 273.15;

% Saturation vapor pressures from XSteam
% XSteam psat_T input: temperature [degC]
% XSteam psat_T output: pressure [bar]
p_sat_H2O_ads = XSteam('psat_T', T_ads_C);
p_sat_H2O_des = XSteam('psat_T', T_des_C);

%% -------------------------
% Adsorber H2O profile
% --------------------------
% The adsorber solver is kept on dry-gas basis.
% H2O is added here as a conserved non-adsorbing vapor for RH calculation.

n_H2O_ads = ...
    yH2O_ads_feed_wet / (1 - yH2O_ads_feed_wet) * n_g_feed;

yH2O_ads_in = yH2O_ads_feed_wet;
pH2O_ads_in = yH2O_ads_in * P_bar;
RH_ads_in = pH2O_ads_in / p_sat_H2O_ads;

n_g_ads_wet_profile = ...
    ads_results.n_g_profile + n_H2O_ads;

yH2O_ads_profile = ...
    n_H2O_ads ./ n_g_ads_wet_profile;

pH2O_ads_profile = ...
    yH2O_ads_profile * P_bar;

RH_ads_profile = ...
    pH2O_ads_profile / p_sat_H2O_ads;

%% -------------------------
% Desorber H2O profile
% --------------------------
% Desorber gas phase is CO2 + H2O.
% H2O adsorption is neglected.

yH2O_des_in = 1 - yCO2_steam_in;
pH2O_des_in = yH2O_des_in * P_bar;
RH_des_in = pH2O_des_in / p_sat_H2O_des;

yH2O_des_profile = ...
    1 - des_results.y_profile;

pH2O_des_profile = ...
    yH2O_des_profile * P_bar;

RH_des_profile = ...
    pH2O_des_profile / p_sat_H2O_des;

%% -------------------------
% Print RH results
% --------------------------

fprintf('\nRelative humidity calculation:\n');

fprintf('Adsorber p_sat,H2O               = %.6f bar at %.2f C\n', ...
    p_sat_H2O_ads, T_ads_C);

fprintf('Desorber p_sat,H2O               = %.6f bar at %.2f C\n', ...
    p_sat_H2O_des, T_des_C);

fprintf('\nAdsorber inlet:\n');
fprintf('yH2O                             = %.6f\n', yH2O_ads_in);
fprintf('pH2O                             = %.6f bar\n', pH2O_ads_in);
fprintf('RH                               = %.2f %%\n', 100 * RH_ads_in);

fprintf('\nAdsorber stage outlet RH profile:\n');
for i = 1:N_ads
    fprintf('Stage %d: yH2O = %.6f, pH2O = %.6f bar, RH = %.2f %%\n', ...
        i, yH2O_ads_profile(i), pH2O_ads_profile(i), ...
        100 * RH_ads_profile(i));
end

fprintf('\nDesorber inlet:\n');
fprintf('yH2O                             = %.6f\n', yH2O_des_in);
fprintf('pH2O                             = %.6f bar\n', pH2O_des_in);
fprintf('RH                               = %.2f %%\n', 100 * RH_des_in);

fprintf('\nDesorber stage outlet RH profile:\n');
for i = 1:N_des
    fprintf('Stage %d: yH2O = %.6f, pH2O = %.6f bar, RH = %.2f %%\n', ...
        i, yH2O_des_profile(i), pH2O_des_profile(i), ...
        100 * RH_des_profile(i));
end

fprintf('\nAverage RH:\n');
fprintf('Adsorber average stage-outlet RH = %.2f %%\n', ...
    100 * mean(RH_ads_profile));
fprintf('Desorber average stage-outlet RH = %.2f %%\n', ...
    100 * mean(RH_des_profile));
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
% TSA HEAT DEMAND
% ==========================

% Desorption heat from loading-dependent heat integration
Q_desorption = des_results.Q_des;

% Sorbent sensible heat between adsorber and desorber
Q_sorbent_sensible_gross = ...
    m_s * cp_sorbent * (T_des - T_ads);

Q_sorbent_sensible_net = ...
    (1 - heat_recovery_fraction) * Q_sorbent_sensible_gross;

% Steam generation heat:
% liquid water at 100 C -> saturated steam at 100 C -> ideal-gas steam at T_des
Q_steam_generation = ...
    m_steam_in * ...
    (h_vap_H2O_100C + cp_H2O_vapor_ideal * (T_des - T_feed_water));

% Total external heat input
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

    q_ads_cooling_specific = ...
        1e-6 * ads_results.Q_ads / ads_results.m_CO2_captured;

else

    q_desorption_specific = NaN;
    q_sorbent_sensible_gross_specific = NaN;
    q_sorbent_sensible_net_specific = NaN;
    q_steam_generation_specific = NaN;
    q_total_heat_specific = NaN;
    q_ads_cooling_specific = NaN;

end

fprintf('\nTSA heat demand:\n');

fprintf('Desorption heat                  = %.2f W\n', Q_desorption);
fprintf('Specific desorption heat         = %.4f MJ/kg CO2\n', ...
    q_desorption_specific);

fprintf('Sorbent sensible heat, gross     = %.2f W\n', ...
    Q_sorbent_sensible_gross);
fprintf('Specific gross sensible heat     = %.4f MJ/kg CO2\n', ...
    q_sorbent_sensible_gross_specific);

fprintf('Lean-rich heat recovery fraction = %.2f %%\n', ...
    100 * heat_recovery_fraction);

fprintf('Sorbent sensible heat, net       = %.2f W\n', ...
    Q_sorbent_sensible_net);
fprintf('Specific net sensible heat       = %.4f MJ/kg CO2\n', ...
    q_sorbent_sensible_net_specific);

fprintf('Steam generation heat            = %.2f W\n', ...
    Q_steam_generation);
fprintf('Specific steam generation heat   = %.4f MJ/kg CO2\n', ...
    q_steam_generation_specific);

fprintf('Total heat demand                = %.2f W\n', Q_total_heat);
fprintf('Specific total heat demand       = %.4f MJ/kg CO2\n', ...
    q_total_heat_specific);

fprintf('\nAdsorber cooling duty, separate:\n');
fprintf('Adsorber cooling duty            = %.2f W\n', ads_results.Q_ads);
fprintf('Specific adsorber cooling duty   = %.4f MJ/kg CO2\n', ...
    q_ads_cooling_specific);

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
