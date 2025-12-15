function [reactions, stoich, x_init] = init_gillespie( ...
    N_E, N_I, w_0, w_tot, h, alpha, act, active_init)
% init_gillespie
%
% N_E, N_I      : number of excitatory / inhibitory units
% w_0, w_tot    : weight parameters
% h             : external drive
% alpha         : decay rate
% act           : activation function handle
% active_init   : fraction initially active (0–1)
%
% Outputs
% reactions : cell array of reaction rate functions
% stoich    : stoichiometry matrix
% x_init    : initial state vector
% phi       : input field function handle

    N = N_E + N_I;

    % ---------- weights ----------
    w_E = (w_0 + w_tot) / 2;
    w_I = w_E - w_0;

    % ---------- field ----------
    phi = @(x) ...
        w_E / N_E * sum(x(2:2:N_E*2)) ...
      - w_I / N_I * sum(x(N_E*2+2:2:end)) ...
      + h;

    % ---------- initialize ----------
    reactions = cell(N * 2, 1);
    stoich    = -eye(N * 2);
    x_init    = zeros(N * 2, 1);

    % ---------- excitatory ----------
    for k = 1:N_E
        kk = k;

        reactions{kk*2-1} = @(x) x(2*kk-1) * act(phi(x));
        reactions{kk*2}   = @(x) x(2*kk)   * alpha;

        stoich(kk*2,   kk*2-1) = 1;
        stoich(kk*2-1, kk*2)   = 1;

        if kk <= active_init * N_E
            x_init(kk*2) = 1;
        else
            x_init(kk*2-1) = 1;
        end
    end

    % ---------- inhibitory ----------
    for k = N_E+1:N
        kk = k;

        reactions{kk*2-1} = @(x) x(2*kk-1) * act(phi(x));
        reactions{kk*2}   = @(x) x(2*kk)   * alpha;

        stoich(kk*2,   kk*2-1) = 1;
        stoich(kk*2-1, kk*2)   = 1;

        if kk - N_E <= active_init * N_I
            x_init(kk*2) = 1;
        else
            x_init(kk*2-1) = 1;
        end
    end
end
