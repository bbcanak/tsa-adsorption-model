function beta_eq = Toth_Lewatit_Sutanto(pCO2_bar, T)
    % Toth isotherm for Lewatit VP OC 1065
    % Source: Sutanto et al., 2017, Separation and Purification Technology 184,
    %         12-25.
    
    % pressure in bar
    % temperature in K
    % output in mol/kg

    R = 8.31447;        % J/mol/K

    % Toth parameters
    ns0   = 3.70;     % mol/kg
    chi   = 0.00;
    Tref  = 353.0;    % K
    b0    = 188.6;    % 1/bar
    dH    = 111.0e3;  % J/mol
    t0    = 0.30;
    alpha = 0.5;
    
    % Temperature-dependent parameters
    b  = b0 * exp((dH/(R*Tref)) * (Tref/T - 1));
    ns = ns0 * exp(chi * (1 - T/Tref));
    t  = t0 + alpha * (1 - Tref/T);

    % Toth isotherm
    beta_eq = (ns * b * pCO2_bar) / ((1 + (b * pCO2_bar)^t)^(1/t));
end