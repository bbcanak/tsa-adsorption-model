function dH_ads = Toth_heat_adsorption(beta, T, Isotherm)
% Toth_heat_adsorption
% Calculates the loading-dependent magnitude of the CO2 heat of adsorption
% from the Toth model.
%
% Inputs:
%   beta      CO2 loading [mol/kg]
%   T         temperature [K]
%   Isotherm  name of selected Toth isotherm
%
% Output:
%   dH_ads    heat of adsorption magnitude [J/mol CO2]

R = 8.31447;    % J/mol/K

%% =========================
% SELECT TOTH PARAMETERS
% ==========================

switch Isotherm

    case 'Toth_Lewatit_Zerobin'

        ns0   = 3.13;       % mol/kg
        chi   = 0.00;
        Tref  = 343.0;      % K
        dH0   = 106.0e3;    % J/mol
        t0    = 0.34;
        alpha = 0.42;

    case 'Toth_Lewatit_Sutanto'

        ns0   = 3.70;       % mol/kg
        chi   = 0.00;
        Tref  = 353.0;      % K
        dH0   = 111.0e3;    % J/mol
        t0    = 0.30;
        alpha = 0.50;

    case 'Toth_Lewatit_Veneman'

        ns0   = 3.40;       % mol/kg
        chi   = 0.00;
        Tref  = 353.0;      % K
        dH0   = 86.70e3;    % J/mol
        t0    = 0.30;
        alpha = 0.14;

    case 'Toth_13X_Zerobin'

        ns0   = 19.0;       % mol/kg
        chi   = 0.00;
        Tref  = 343.0;      % K
        dH0   = 61.7e3;     % J/mol
        t0    = 0.30;
        alpha = 0.49;

    otherwise

        error('Unknown isotherm model: %s', Isotherm);

end

%% =========================
% TEMPERATURE-DEPENDENT TOTH PARAMETERS
% ==========================

ns = ns0 * exp(chi * (1 - T / Tref));

t = t0 + alpha * (1 - Tref / T);

%% =========================
% NORMALIZED LOADING
% ==========================

theta = beta / ns;

% Numerical protection
theta = max(theta, 1e-12);
theta = min(theta, 1 - 1e-12);

%% =========================
% TOTH HEAT OF ADSORPTION
% ==========================

theta_t = theta.^t;

term_1 = log(1 - theta_t) / t^2;

term_2 = theta_t .* log(theta) ./ ...
    (t .* (1 - theta_t));

dH_ads = dH0 + ...
    R * alpha * Tref .* (term_1 + term_2);

end