function [reactions, stoich, x_init] = init_simplified( ...
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
    phi = @(x) w_E / N_E * x(1) - w_I / N_I * x(2) + h;

    % ---------- initialize ----------
    reactions = {@(x) (N_E - x(1)) * act(phi(x)), ...
        @(x) x(1) * alpha, ...
        @(x) (N_I - x(2)) * act(phi(x)), ...
        @(x) x(2) * alpha};
    stoich    = [1 -1 0 0; 0 0 1 -1];
    x_init    = [round(N_E * active_init); round(N_I * active_init)];

end
