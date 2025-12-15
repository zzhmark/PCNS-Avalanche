%% Definition
% Here we run the simplified model for figure 3 data and define the
% parameters as the paper suggests.

% from fig3
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

% fig7
conn_pct = 0.17;    % sparsity
wff = 10;
fname = 'fig7/data.mat';

if ~exist('fig7','dir'); mkdir(path); end

%% Init
[MD,~] = qr(randn(N));   % orthogonal
[MS,~] = qr(randn(N));   % orthogonal

k = round(conn_pct * N);

CD = zeros(N);
CS = zeros(N);

idx = randperm(N,k);

% here we simply use identical values for the eigenvalues, which are
% comparable to previous experiments
CD(idx,idx) = w_0;
CS(idx,idx) = wff;

Wdiff = MD' * CD * MD;   % W_E - W_I
Wsum  = MS' * CS * MS;   % W_E + W_I

WE = (Wsum + Wdiff)/2;
WI = (Wsum - Wdiff)/2;

% enforce dale's principle
for j = 1:N
    WE(:,j) = max(WE(:,j), 0);   % excitatory
    WI(:,j) = max(WI(:,j), 0);   % inhibitory
end

W = [ WE  -WI; WE  -WI ];
[~, ~, x0] = init_gillespie(N_E, N_I, w_0, wff, h, alpha, act, active_init);

%% Gillespie

rng('shuffle');

t = 0;
x = x0(:);

Tcell = {t};
Xcell = {x.'};

k = 1;
while t < ttot

    % input
    s = W * x + h;

    % propensities
    a = zeros(N,1);
    a(x==0) = act(s(x==0));   % Q -> A
    a(x==1) = alpha;       % A -> Q

    a0 = sum(a);
    if a0 == 0, break; end

    % Gillespie step
    tau = -log(rand()) / a0;
    r   = rand() * a0;
    i   = find(cumsum(a) >= r, 1);

    % update state
    x(i) = 1 - x(i);
    t = t + tau;

    % store
    k = k + 1;
    Tcell{k,1} = t;
    Xcell{k,1} = x.';
    msg = sprintf('\rTime Progress (ms): %f / %f', t, ttot); 
    fprintf([repmat('\b',1,100) msg]);
end
fprintf('\n');

% ---- concatenate once ----
T = vertcat(Tcell{:});   % (N,1)
X = vertcat(Xcell{:});   % (N,numSpecies)

save(fname, 'T', 'X');

%% load saved
load(fname, 'T', 'X');
active = X(:, 2:2:2*N);
D = diff(active,1,1);        % diff over T (rows)
D(D < 0) = 0;           % remove negatives
event = [zeros(1,N); D];        % prepend first row

%% Plots Fig7A

% comptue mean fr
[Tb, nspk] = timebin_sum(T, event, tbin);
fr = nspk / tbin * 1000;    % in Hz
fr = gauss_smooth_1d(fr, sigma / tbin);
mfr = mean(fr, 2);

f = plot_raster_with_mfr(Tb, nspk, mfr, ttot);

exportgraphics(f,"fig7/A.png", "Resolution", 300);

%% Plot Fig7B

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

exportgraphics(f,"fig7/B.png", "Resolution", 300)

%% Plot Fig7C
w_E = (w_0 + wff) / 2;
w_I = w_E - w_0;
f = plot_phase_plane(X, T, tbin, alpha, w_E, w_I, h, act, N_E, N_I);
exportgraphics(f,"fig7/C.png", "Resolution", 300);