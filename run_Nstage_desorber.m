%% run_Nstage_desorber.m
% General N-stage counter-current equilibrium desorber model
% CO2 desorption from loaded sorbent using selectable Toth isotherm

clear; clc;

%% =========================
% USER INPUTS
% ==========================

% Number of desorber stages
N = 4;

% Sorbent / isotherm model
Isotherm = 'Toth_Lewatit_Zerobin';
%Isotherm = 'Toth_Lewatit_Veneman';
%Isotherm = 'Toth_Lewatit_Sutanto';
%Isotherm = 'Toth_13X_Zerobin';


% Operating conditions
T_des = 120 + 273.15;    % K
P_bar = 1.01325;         % bar

% Solid inlet to desorber
m_s = 0.10;              % kg/s sorbent circulation rate
beta_rich = 1.2;        % mol CO2/kg sorbent entering desorber

% Stripping steam inlet
n_steam_in = 0.20;       % mol H2O/s
yCO2_steam_in = 0.0;     % pure stripping steam

%% =========================
% SOLVE N-STAGE DESORBER
% ==========================

results = solve_desorber_Nstage(N, n_steam_in, yCO2_steam_in, ...
                                m_s, beta_rich, ...
                                T_des, P_bar, ...
                                Isotherm);

%% =========================
% DISPLAY RESULTS
% ==========================

fprintf('\n===== N-STAGE COUNTER-CURRENT DESORBER RESULTS =====\n');
fprintf('Number of stages          = %d\n', N);
fprintf('Isotherm                  = %s\n', Isotherm);
fprintf('T_des                     = %.2f K\n', T_des);
fprintf('P                         = %.5f bar\n', P_bar);
fprintf('Stripping steam flow      = %.6f mol H2O/s\n', n_steam_in);
fprintf('Steam inlet yCO2          = %.6f\n', yCO2_steam_in);
fprintf('Sorbent flow              = %.6f kg/s\n', m_s);
fprintf('Rich loading              = %.6f mol/kg\n', beta_rich);

fprintf('\nOverall results:\n');
fprintf('CO2-rich gas yCO2 out     = %.6f\n', results.y_CO2_rich_gas_out);
fprintf('Lean loading              = %.6f mol/kg\n', results.beta_lean);
fprintf('Working capacity          = %.6f mol/kg\n', results.working_capacity);
fprintf('CO2 desorbed              = %.6f mol/s\n', results.n_CO2_desorbed);
fprintf('CO2 desorbed              = %.6f kg/s\n', results.m_CO2_desorbed);
fprintf('Desorber heat duty        = %.2f W\n', results.Q_des);
fprintf('Specific desorption duty  = %.4f MJ/kg CO2\n', results.q_des_specific);

fprintf('\nStagewise gas CO2 mole fractions:\n');
disp(results.y_profile.');

fprintf('Stagewise sorbent loadings, mol/kg:\n');
disp(results.beta_profile.');

fprintf('\nStagewise desorber heat duties:\n');
for i = 1:N
    fprintf('Stage %d: CO2 desorbed = %.8f mol/s, Q_des = %.2f W\n', ...
        i, results.n_CO2_desorbed_stage(i), results.Q_des_stage(i));
end
%% =========================
% DESORBER MASS BALANCE CHECK
% ==========================

CO2_by_solid_des = m_s * ...
    (results.beta_rich - results.beta_lean);

fprintf('\nDesorber mass balance check:\n');

fprintf('CO2 desorbed to gas             = %.8f mol/s\n', ...
    results.n_CO2_desorbed);

fprintf('CO2 released by solid           = %.8f mol/s\n', ...
    CO2_by_solid_des);

fprintf('Gas-solid mismatch              = %.8e mol/s\n', ...
    results.n_CO2_desorbed - CO2_by_solid_des);
%% =========================
% PLOTS
% ==========================

stage = 1:N;

figure;
plot(stage, results.y_profile, 'o-', 'LineWidth', 1.5);
xlabel('Stage number');
ylabel('Gas-phase y_{CO2}');
title('Stagewise Gas CO2 Profile in Desorber');
grid on;

figure;
plot(stage, results.beta_profile, 's-', 'LineWidth', 1.5);
xlabel('Stage number');
ylabel('\beta_{CO2} [mol/kg]');
title('Stagewise Sorbent Loading Profile in Desorber');
grid on;

figure;
bar(stage, results.Q_des_stage/1000);
xlabel('Stage number');
ylabel('Desorber heat duty [kW]');
title('Stagewise Desorber Heat Duty');
grid on;
