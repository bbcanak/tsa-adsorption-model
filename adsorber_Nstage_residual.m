function res = adsorber_Nstage_residual(x, N, ...
                                         n_g_feed, yCO2_feed, ...
                                         m_s, beta_lean, ...
                                         T_ads, P_bar, ...
                                         Isotherm)
% adsorber_Nstage_residual
% Residual equations for N-stage counter-current equilibrium adsorber.
%
% Unknowns:
%   y(i)     gas CO2 mole fraction leaving stage i
%   beta(i)  solid CO2 loading leaving stage i
%
% Stage indexing:
%   Gas flows from stage 1 to stage N.
%   Solid flows from stage N to stage 1.
%
% Therefore:
%   Gas inlet to stage 1 = feed gas
%   Gas inlet to stage i = gas outlet from stage i-1
%
%   Solid inlet to stage N = lean sorbent
%   Solid inlet to stage i = solid outlet from stage i+1

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
    % STAGEWISE EQUATIONS
    % ==========================

    n_g_out_previous = NaN;

    for i = 1:N

        %% -------------------------
        % Gas inlet to stage i
        % --------------------------

        if i == 1
            n_g_in_i = n_g_feed;
            y_in_i = yCO2_feed;
        else
            y_in_i = y(i-1);
            n_g_in_i = n_g_out_previous;
        end

        %% -------------------------
        % Gas outlet from stage i
        % --------------------------

        y_out_i = y(i);

        % Inert balance:
        n_g_out_i = n_g_in_i * (1 - y_in_i) / (1 - y_out_i);

        %% -------------------------
        % Solid inlet to stage i
        % --------------------------

        if i == N
            beta_in_i = beta_lean;
        else
            beta_in_i = beta(i+1);
        end

        %% -------------------------
        % Solid outlet from stage i
        % --------------------------

        beta_out_i = beta(i);

        %% -------------------------
        % CO2 mass balance
        % --------------------------

        res(i) = n_g_in_i * y_in_i + m_s * beta_in_i ...
               - n_g_out_i * y_out_i - m_s * beta_out_i;

        %% -------------------------
        % Equilibrium condition
        % --------------------------

        pCO2_out_i = y_out_i * P_bar;
        beta_eq_i = feval(Isotherm, pCO2_out_i, T_ads);

        res(N+i) = beta_out_i - beta_eq_i;

        %% -------------------------
        % Store gas outlet flow
        % --------------------------

        n_g_out_previous = n_g_out_i;
    end
end