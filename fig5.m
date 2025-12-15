%% Definition
% Here we run the simplified model for figure 5 data and define the
% parameters as the paper suggests.

% params from fig3
h = 0.001;     % external input for excitatory and inhibitory
alpha = 0.1;    % propensity for active to inactive, ms, from this we also assume that the max FR of a neuron is 100Hz
w_0 = 0.2;       % the difference in weights between excitatory and inhibitory
r_E = 0.5;       % we assume that the excitatory neurons makes up 50%
sigma = 5;      % gaussian sigma for smoothing in ms
act = @(s) (s > 0) .* tanh(s);  % activation function
active_init = 0.5;  % proportion of active neurons in initiailization
ttot = 1000;    % total simulation time in ms
x_min = 10;     % power law fit parameter
n_ISI = 50;         % # neuron sampled for ISI distribution
tbin = 0.1;     % time bin in ms

% params for multiple wff simulation
wff = 10;

if ~exist('fig5','dir'); mkdir(path); end
if ~exist('fig5_simplified','dir'); mkdir(path); end

%% Gillespie (very slow, and the output data is very big, so saving X as binary)

sizes = [2e3 5e3];        % total population, only able to do these due to limited memory and time
for N = sizes
    fname = fullfile('fig5', [num2str(N) '.mat']);
    N_E = round(N * r_E);  % the exact # of excitatory neurons
    N_I = N - N_E; % the exact # of excitatory neurons
    [reactions, stoich, x_init] = init_gillespie(N_E, N_I, w_0, wff, h, alpha, act, active_init);
    [T, X] = gillespie_binary(reactions, stoich, x_init, ttot);
    save(fname, 'T', 'X', '-v7.3');     % X can be very large, so use a different saving version
end

%% Plot fig5 raster plots
files = dir(fullfile("fig5", "*.mat"));   % list .mat files in folder 

for k = 1:numel(files)
    fname = files(k).name;                 % e.g. '36.mat'
    name_no_ext = erase(fname, ".mat");    % strip extension
    N = str2double(name_no_ext);         % convert to number
    load(fullfile('fig5', fname));     % load variables 
    active = X(:, 2:2:end);
    D = diff(active,1,1);        % diff over T (rows)
    D(D < 0) = 0;           % remove negatives
    event = [zeros(1,N); D];        % prepend first row
    
    % comptue mean fr
    [Tb, nspk] = timebin_sum(T, event, tbin);
    fr = nspk / tbin * 1000;    % in Hz
    fr = gauss_smooth_1d(fr, sigma / tbin);
    mfr = mean(fr, 2);
    
    f = plot_raster_with_mfr(Tb, nspk, mfr, ttot);
    
    fname = fullfile('fig5', [num2str(N) '.png']);
    exportgraphics(f, fname, "Resolution", 300);
end

%% Try simplified model
sizes = [2e3 5e3 1e4 1e5];        % total population
for N = sizes
    fname = fullfile('fig5_simplified', [num2str(N) '.mat']);
    N_E = round(N * r_E);  % the exact # of excitatory neurons
    N_I = N - N_E; % the exact # of excitatory neurons
    [reactions, stoich, x_init] = init_simplified(N_E, N_I, w_0, wff, h, alpha, act, active_init);
    [T, X] = gillespie(reactions, stoich, x_init, ttot);
    save(fname, 'T', 'X');
end

%% Plot fig5 raster plots
files = dir(fullfile("fig5_simplified", "*.mat"));   % list .mat files in folder 
files = natsortfiles(files);
f = figure;
hold on
for k = 1:numel(files)
    fname = files(k).name;                 % e.g. '36.mat'
    name_no_ext = erase(fname, ".mat");    % strip extension
    N = str2double(name_no_ext);         % convert to number
    load(fullfile('fig5_simplified', fname));     % load variables 
    active = sum(X, 2) / N;
    plot(T, active + 4 - k, 'DisplayName', name_no_ext);
end
legend show 
xlabel('Time (ms)');
ylabel('Firing proportions');
yticks([])
yticklabels([])
xlim([0 ttot]);
hold off
fname = fullfile('fig5_simplified', 'combined.png');
exportgraphics(f, fname, "Resolution", 300);