%% run_Nstage_adsorber.m
% General N-stage counter-current equilibrium adsorber model
% CO2 adsorption from flue gas using selectable Toth isotherm

clear; clc;

%% =========================
% USER INPUTS
% ==========================

% Number of stages
N = 5;        % First test N = 2, then use N = 5

% Sorbent / isotherm model
Isotherm = 'Toth_Lewatit_Zerobin';
%Isotherm = 'Toth_Lewatit_Veneman';
%Isotherm = 'Toth_Lewatit_Sutanto';
%Isotherm = 'Toth_13X_Zerobin';


% Operating conditions
T_ads = 50 + 273.15;     % K
P_bar = 1.01325;         % bar

% Gas inlet to adsorber
n_g_feed = 1.0;          % mol/s total gas feed
yCO2_feed = 0.1;        % mol fraction CO2 in feed gas

% Solid inlet to adsorber
m_s = 0.1;              % kg/s sorbent circulation rate
beta_lean = 0.00200828;        % mol CO2/kg sorbent entering adsorber

%% =========================
% SOLVE N-STAGE ADSORBER
% ==========================

results = solve_adsorber_Nstage(N, n_g_feed, yCO2_feed, ...
    m_s, beta_lean, ...
    T_ads, P_bar, ...
    Isotherm);

%% =========================
% DISPLAY RESULTS
% ==========================

fprintf('\n===== N-STAGE COUNTER-CURRENT ADSORBER RESULTS =====\n');
fprintf('Number of stages          = %d\n', N);
fprintf('Isotherm                  = %s\n', Isotherm);
fprintf('T_ads                     = %.2f K\n', T_ads);
fprintf('P                         = %.5f bar\n', P_bar);
fprintf('Gas feed flow             = %.6f mol/s\n', n_g_feed);
fprintf('Feed yCO2                 = %.6f\n', yCO2_feed);
fprintf('Sorbent flow              = %.6f kg/s\n', m_s);
fprintf('Lean loading              = %.6f mol/kg\n', beta_lean);

fprintf('\nOverall results:\n');
fprintf('Outlet yCO2               = %.6f\n', results.y_gas_out);
fprintf('Rich loading              = %.6f mol/kg\n', results.beta_rich);
fprintf('CO2 captured              = %.6f mol/s\n', results.n_CO2_captured);
fprintf('CO2 captured              = %.6f kg/s\n', results.m_CO2_captured);
fprintf('Capture fraction          = %.2f %%\n', 100 * results.capture_fraction);
fprintf('Adsorber cooling duty     = %.2f W\n', results.Q_ads);
fprintf('Specific cooling duty     = %.4f MJ/kg CO2\n', results.q_ads_specific);

fprintf('\nStagewise gas CO2 mole fractions:\n');
disp(results.y_profile.');

fprintf('Stagewise sorbent loadings, mol/kg:\n');
disp(results.beta_profile.');
fprintf('\nStagewise adsorber cooling duties:\n');
for i = 1:N
    fprintf('Stage %d: CO2 captured = %.8f mol/s, Q_ads = %.2f W\n', ...
        i, results.n_CO2_captured_stage(i), results.Q_ads_stage(i));
end
%% =========================
% ADSORBER MASS BALANCE CHECK
% ==========================

CO2_by_solid_ads = m_s * ...
    (results.beta_rich - results.beta_lean);

fprintf('\nAdsorber mass balance check:\n');

fprintf('CO2 captured from gas            = %.8f mol/s\n', ...
        results.n_CO2_captured);

fprintf('CO2 uptake by solid              = %.8f mol/s\n', ...
        CO2_by_solid_ads);

fprintf('Gas-solid mismatch               = %.8e mol/s\n', ...
        results.n_CO2_captured - CO2_by_solid_ads);
%% =========================
% PLOTS
% ==========================

stage = 1:N;

figure;
plot(stage, results.y_profile, 'o-', 'LineWidth', 1.5);
xlabel('Stage number');
ylabel('Gas-phase y_{CO2}');
title('Stagewise Gas CO2 Profile in Adsorber');
grid on;

figure;
plot(stage, results.beta_profile, 's-', 'LineWidth', 1.5);
xlabel('Stage number');
ylabel('\beta_{CO2} [mol/kg]');
title('Stagewise Sorbent Loading Profile in Adsorber');
grid on;

figure;
bar(stage, results.Q_ads_stage);
xlabel('Stage number');
ylabel('Adsorber cooling duty [W]');
title('Stagewise Adsorber Cooling Duty');
grid on;
