%% Definition
% Here we run the simplified model for figure 3 data and define the
% parameters as the paper suggests.

N = 800;        % total population
h = 0.001;     % external input for excitatory and inhibitory
w_0 = 0.2;       % the difference in weights between excitatory and inhibitory
alpha = 0.1;    % propensity for active to inactive, ms, from this we also assume that the max FR of a neuron is 100Hz
r_E = 0.5;       % we assume that the excitatory neurons makes up 50%
tbin = 0.1;     % time bin in ms
sigma = 5;      % gaussian sigma for smoothing in ms
act = @(s) (s > 0) .* tanh(s);  % activation function
active_init = 0.5;  % proportion of active neurons in initiailization
N_E = round(N * r_E);  % the exact # of excitatory neurons
N_I = N - N_E; % the exact # of excitatory neurons
ttot = 1000;    % total simulation time in ms
x_min = 10;     % power law fit parameter
n_ISI = 50;         % # neuron sampled for ISI distribution

if ~exist('fig3','dir'); mkdir(path); end

%% ========================================= Simulation for wE + wI = 0.8
w_tot = 0.8;
fname = 'fig3/data_08.mat';
w_E = (w_0 + w_tot) / 2;
w_I = w_E - w_0;
[reactions, stoich, x_init] = init_gillespie(N_E, N_I, w_0, w_tot, h, alpha, act, active_init);

% run Gillespie (skip to load data)
if isfile(fname)
    choice = questdlg(sprintf('Result "%s" exists. Overwrite or Load?', fname), ...
                       'Confirm Overwrite', 'Overwrite','Load','Cancel','Load');
else
    choice = 'Overwrite';
end
switch choice
    case 'Overwrite'
        [T, X] = gillespie(reactions, stoich, x_init, ttot);
        save(fname, 'T', 'X');
    case 'Load'
        load(fname);
end

active = X(:, 2:2:2*N);
D = diff(active,1,1);        % diff over T (rows)
D(D < 0) = 0;           % remove negatives
event = [zeros(1,N); D];        % prepend first row

%% Plot Fig3A (load above data before run)

% comptue mean fr
[Tb, nspk] = timebin_sum(T, event, tbin);
fr = nspk / tbin * 1000;    % in Hz
fr = gauss_smooth_1d(fr, sigma / tbin);
mfr = mean(fr, 2);

f = plot_raster_with_mfr(Tb, nspk, mfr, ttot);

exportgraphics(f,"fig3/A.png", "Resolution", 300);

%% Plot Fig3D (load above data before run)

active = X(:, 2:2:2*N);
D = diff(active,1,1);        % diff over T (rows)
D(D < 0) = 0;           % remove negatives
event = [zeros(1,N); D];        % prepend first row

% compute ISI (network level)
isi_network = compute_isi(T, max(event, [], 2));
dt = mean(isi_network);

% compute avalanche sizes
avalanche_sizes = avalanche_sizes_from_isi(isi_network, dt);

% geometric fitting, MLE
p_hat = 1 / (1 + mean(avalanche_sizes));
mu = -log(p_hat) / dt;

% power law fitting, MLE
avalanche_sizes2 = avalanche_sizes(avalanche_sizes > x_min);
beta = 1 + numel(avalanche_sizes2) / (sum(log(avalanche_sizes2 / (x_min - 0.5))));

% exponential fitting for ISI
rng(41);
isi_neurons = cell(n_ISI);
p = randperm(N, n_ISI);
for i = 1:n_ISI
    isi_neurons{i} = compute_isi(T, event(:, p(i)));
end
isi_neurons = vertcat(isi_neurons{:});
lambda_hat = 1 / mean(isi_neurons);

f = plot_avalanche_isi_fit(avalanche_sizes, avalanche_sizes2, isi_neurons, p_hat, beta, lambda_hat, dt);

exportgraphics(f,"fig3/D.png", "Resolution", 300)
%% Plot Fig3G (load above data before run)

f = plot_phase_plane(X, T, tbin, alpha, w_E, w_I, h, act, N_E, N_I);
exportgraphics(f,"fig3/G.png", "Resolution", 300)

%% ========================================= Simulation for wE + wI = 1.8
w_tot = 1.8;
fname = 'fig3/data_18.mat';
w_E = (w_0 + w_tot) / 2;
w_I = w_E - w_0;
[reactions, stoich, x_init] = init_gillespie(N_E, N_I, w_0, w_tot, h, alpha, act, active_init);

% run Gillespie (skip to load data)
if isfile(fname)
    choice = questdlg(sprintf('Result "%s" exists. Overwrite or Load?', fname), ...
                       'Confirm Overwrite', 'Overwrite','Load','Cancel','Load');
else
    choice = 'Overwrite';
end
switch choice
    case 'Overwrite'
        [T, X] = gillespie(reactions, stoich, x_init, ttot);
        save(fname, 'T', 'X');
    case 'Load'
        load(fname);
end

active = X(:, 2:2:2*N);
D = diff(active,1,1);        % diff over T (rows)
D(D < 0) = 0;           % remove negatives
event = [zeros(1,N); D];        % prepend first row

%% Plot Fig3B (load above data before run)

% comptue mean fr
[Tb, nspk] = timebin_sum(T, event, tbin);
fr = nspk / tbin * 1000;
fr = gauss_smooth_1d(fr, sigma / tbin);
mfr = mean(fr, 2);

f = plot_raster_with_mfr(Tb, nspk, mfr, ttot);

exportgraphics(f,"fig3/B.png", "Resolution", 300);

%% Plot Fig3E (load above data before run)

% compute ISI (network level)
isi_network = compute_isi(T, max(event, [], 2));
dt = mean(isi_network);

% compute avalanche sizes
avalanche_sizes = avalanche_sizes_from_isi(isi_network, dt);

% geometric fitting, MLE
p_hat = 1 / (1 + mean(avalanche_sizes));
mu = -log(p_hat) / dt;

% power law fitting, MLE
avalanche_sizes2 = avalanche_sizes(avalanche_sizes > x_min);
beta = 1 + numel(avalanche_sizes2) / (sum(log(avalanche_sizes2 / (x_min - 0.5))));

% exponential fitting for ISI
rng(41);
isi_neurons = cell(n_ISI);
p = randperm(N, n_ISI);
for i = 1:n_ISI
    isi_neurons{i} = compute_isi(T, event(:, p(i)));
end
isi_neurons = vertcat(isi_neurons{:});
lambda_hat = 1 / mean(isi_neurons);

f = plot_avalanche_isi_fit(avalanche_sizes, avalanche_sizes2, isi_neurons, p_hat, beta, lambda_hat, dt);

exportgraphics(f,"fig3/E.png", "Resolution", 300)
%% Plot Fig3H (load above data before run)

f = plot_phase_plane(X, T, tbin, alpha, w_E, w_I, h, act, N_E, N_I);
exportgraphics(f,"fig3/H.png", "Resolution", 300)

%% ======================================== Simulation for wE + wI = 13.8

w_tot = 13.8;
fname = 'fig3/data_138.mat';
w_E = (w_0 + w_tot) / 2;
w_I = w_E - w_0;
[reactions, stoich, x_init] = init_gillespie(N_E, N_I, w_0, w_tot, h, alpha, act, active_init);

% run Gillespie (skip to load data)
if isfile(fname)
    choice = questdlg(sprintf('Result "%s" exists. Overwrite or Load?', fname), ...
                       'Confirm Overwrite', 'Overwrite','Load','Cancel','Load');
else
    choice = 'Overwrite';
end
switch choice
    case 'Overwrite'
        [T, X] = gillespie(reactions, stoich, x_init, ttot);
        save(fname, 'T', 'X');
    case 'Load'
        load(fname);
end

active = X(:, 2:2:2*N);
D = diff(active,1,1);        % diff over T (rows)
D(D < 0) = 0;           % remove negatives
event = [zeros(1,N); D];        % prepend first row


%% Plot Fig3C (load above data before run)

% comptue mean fr
[Tb, nspk] = timebin_sum(T, event, tbin);
fr = nspk / tbin * 1000;
fr = gauss_smooth_1d(fr, sigma / tbin);
mfr = mean(fr, 2);

f = plot_raster_with_mfr(Tb, nspk, mfr, ttot);

exportgraphics(f,"fig3/C.png", "Resolution", 300)

%% Plot Fig3F (load above data before run)

% compute ISI (network level)
isi_network = compute_isi(T, max(event, [], 2));
dt = mean(isi_network);

% compute avalanche sizes
avalanche_sizes = avalanche_sizes_from_isi(isi_network, dt);

% geometric fitting, MLE
p_hat = 1 / (1 + mean(avalanche_sizes));
mu = -log(p_hat) / dt;

% power law fitting, MLE
avalanche_sizes2 = avalanche_sizes(avalanche_sizes > x_min);
beta = 1 + numel(avalanche_sizes2) / (sum(log(avalanche_sizes2 / (x_min - 0.5))));

% exponential fitting for ISI
rng(41);
isi_neurons = cell(n_ISI);
p = randperm(N, n_ISI);
for i = 1:n_ISI
    isi_neurons{i} = compute_isi(T, event(:, p(i)));
end
isi_neurons = vertcat(isi_neurons{:});
lambda_hat = 1 / mean(isi_neurons);

f = plot_avalanche_isi_fit(avalanche_sizes, avalanche_sizes2, isi_neurons, p_hat, beta, lambda_hat, dt);

exportgraphics(f,"fig3/F.png", "Resolution", 300)

%% Plot Fig3I (load above data before run)

f = plot_phase_plane(X, T, tbin, alpha, w_E, w_I, h, act, N_E, N_I);
exportgraphics(f,"fig3/I.png", "Resolution", 300)