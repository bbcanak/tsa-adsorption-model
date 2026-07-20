function res = desorber_Nstage_residual(x, N, ...
                                         n_steam_in, yCO2_steam_in, ...
                                         m_s, beta_rich, ...
                                         T_des, P_bar, ...
                                         Isotherm)
% desorber_Nstage_residual
% Residual equations for N-stage counter-current equilibrium desorber.
%
% Unknowns:
%   y(i)     gas CO2 mole fraction leaving stage i
%   beta(i)  solid CO2 loading leaving stage i
%
% Stage indexing:
%   Solid flows from Stage 1 to Stage N.
%   Gas flows from Stage N to Stage 1.
%
% Therefore:
%   Solid inlet to Stage 1 = beta_rich
%   Solid inlet to Stage i = beta outlet from Stage i-1
%
%   Stripping steam enters Stage N
%   Gas inlet to Stage i = Gas outlet from Stage i+1

    %% =========================
    % UNPACK VARIABLES
    % ==========================

    y = x(1:N);
    beta = x(N+1:2*N);

    res = zeros(2*N, 1);

   %% =========================
% PHYSICAL DOMAIN CHECK
% ==========================

if any(~isreal(y)) || any(~isreal(beta)) || ...
   any(y <= 0) || any(y >= 1) || any(beta < 0)

    res(:) = 1e3;
    return;
end
    %% =========================
    % GAS FLOW CALCULATION
    % ==========================
    % Gas flows from Stage N to Stage 1.
    % So calculate gas outlet flows backwards first.

    n_g_out = zeros(1, N);

    for i = N:-1:1

        % Gas inlet to stage i
        if i == N
            n_g_in_i = n_steam_in;
            y_in_i = yCO2_steam_in;
        else
            n_g_in_i = n_g_out(i+1);
            y_in_i = y(i+1);
        end

        % Gas outlet from stage i
        y_out_i = y(i);

        % H2O component balance
        % H2O adsorption is neglected; therefore the steam molar flow
        % is conserved across each desorber stage.
        n_g_out_i = n_g_in_i * (1 - y_in_i) / (1 - y_out_i);

        n_g_out(i) = n_g_out_i;
    end

    %% =========================
    % STAGEWISE EQUATIONS
    % ==========================

    for i = 1:N

        %% -------------------------
        % Solid inlet to stage i
        % --------------------------

        if i == 1
            beta_in_i = beta_rich;
        else
            beta_in_i = beta(i-1);
        end

        beta_out_i = beta(i);

        %% -------------------------
        % Gas inlet to stage i
        % --------------------------

        if i == N
            n_g_in_i = n_steam_in;
            y_in_i = yCO2_steam_in;
        else
            n_g_in_i = n_g_out(i+1);
            y_in_i = y(i+1);
        end

        %% -------------------------
        % Gas outlet from stage i
        % --------------------------

        n_g_out_i = n_g_out(i);
        y_out_i = y(i);

        %% -------------------------
        % CO2 mass balance
        % --------------------------
        % CO2 transfers from solid to gas:
        % gas in + solid in = gas out + solid out

        res(i) = n_g_in_i * y_in_i + m_s * beta_in_i ...
               - n_g_out_i * y_out_i - m_s * beta_out_i;

        %% -------------------------
        % Equilibrium condition
        % --------------------------

        pCO2_out_i = y_out_i * P_bar;
        beta_eq_i = feval(Isotherm, pCO2_out_i, T_des);

        res(N+i) = beta_out_i - beta_eq_i;
    end
end
