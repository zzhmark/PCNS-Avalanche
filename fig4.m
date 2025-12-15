%% Definition
% Here we run the simplified model for figure 4 data and define the
% parameters as the paper suggests.

% params from fig3
N = 800;        % total population
h = 0.001;     % external input for excitatory and inhibitory
w_0 = 0.2;       % the difference in weights between excitatory and inhibitory
alpha = 0.1;    % propensity for active to inactive, ms, from this we also assume that the max FR of a neuron is 100Hz
r_E = 0.5;       % we assume that the excitatory neurons makes up 50%
sigma = 5;      % gaussian sigma for smoothing in ms
act = @(s) (s > 0) .* tanh(s);  % activation function
active_init = 0.5;  % proportion of active neurons in initiailization
N_E = round(N * r_E);  % the exact # of excitatory neurons
N_I = N - N_E; % the exact # of excitatory neurons
ttot = 1000;    % total simulation time in ms
x_min = 10;     % power law fit parameter
n_ISI = 50;         % # neuron sampled for ISI distribution

% params for multiple wff simulation
wff_step = 0.1;
wff_start = 0.8;
wff_end = 6.4;
tbin = 1;     % time bin in ms

if ~exist('fig4','dir'); mkdir('fig4'); end

%% Gillespie in parallel!
wff = wff_start:wff_step:wff_end;
parfor i=1:numel(wff)
    fname = fullfile('fig4', [num2str(wff(i)) '.mat']);
    [reactions, stoich, x_init] = init_gillespie(N_E, N_I, w_0, wff(i), h, alpha, act, active_init);
    [T, X] = gillespie(reactions, stoich, x_init, ttot);
    s = struct('T', T, 'X', X);
    save(fname, "-fromstruct", s);
end

%% Loading results and compute firing rate
files = dir(fullfile("fig4", "*.mat"));   % list .mat files in folder 

wff = nan(numel(files),1);   % preallocate numeric array

mfr = zeros(numel(files), 1);
sfr = zeros(numel(files), 1);
cv = zeros(numel(files), 1);
for k = 1:numel(files)
    fname = files(k).name;                 % e.g. '36.mat'
    name_no_ext = erase(fname, ".mat");    % strip extension
    num = str2double(name_no_ext);         % convert to number
    wff(k) = num;                         % store result
    
    % optional: load each file
    load(fullfile('fig4', fname));     % load variables 

    % active = X(:, 2:2:2*N_E);
    % D = diff(active,1,1);        % diff over T (rows)
    % D(D < 0) = 0;           % remove negatives
    % spike_E = [zeros(1,N_E); D];        % prepend first row
    % [Tb, nspk_E] = timebin_sum(T, spike_E, tbin);
    % nspk_E = sum(nspk_E, 2);
    % 
    % active = X(:, 2*N_E+2:2:end);
    % D = diff(active,1,1);        % diff over T (rows)
    % D(D < 0) = 0;           % remove negatives
    % spike_I = [zeros(1,N_I); D];        % prepend first row
    % [Tb, nspk_I] = timebin_sum(T, spike_I, tbin);
    % nspk_I = sum(nspk_I, 2);
    % 
    % kt = nspk_I + nspk_E;
    % cv(k) = std(kt) / mean(kt);
    % 
    % % filtered E/I
    % E = zeros(size(nspk_E)); 
    % I = zeros(size(nspk_I));
    % E(1) = 0.1; 
    % I(1) = 0;
    % 
    % for t = 2:numel(nspk_E)
    %     E(t) = (1-alpha*tbin)*E(t-1) + nspk_E(t)/N_E;
    %     I(t) = (1-alpha*tbin)*I(t-1) + nspk_I(t)/N_I;
    % end
    % 
    % fr = (2-E-I).*act(w_E.*E - w_I.*I + h);
    % fr = gauss_smooth_1d(fr, sigma / tbin);

    active = X(:, 2:2:2*N);
    D = diff(active,1,1);        % diff over T (rows)
    D(D < 0) = 0;           % remove negatives
    event = [zeros(1,N); D];        % prepend first row
    [Tb, nspk] = timebin_sum(T, event, tbin);
    kt = mean(nspk, 2);

    cv(k) = std(kt) / mean(kt);
    mfr(k) = mean(kt) / tbin * 1000;
    sfr(k) = std(kt) / tbin * 1000;
end

%% Plot Figure 4A, mean & std of firing rate

f=figure;
[~, idx] = sort(wff);
plot(wff(idx), mfr(idx), 'b-', wff(idx), sfr(idx), 'r-');
legend('mean', 'std')
ylabel('firing rate (Hz/neuron)');
xlabel('w+');

exportgraphics(f,"fig4/A.png", "Resolution", 300)

%% Plot Figure 4B, coefficient of variation

f=figure;
[~, idx] = sort(wff);
plot(wff(idx), cv(idx), 'b-');
ylabel('coefficient of variation');
xlabel('w+');

exportgraphics(f,"fig4/B.png", "Resolution", 300)
