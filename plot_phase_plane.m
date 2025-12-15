function f = plot_phase_plane(X, T, tbin, alpha, w_E, w_I, h, act, N_E, N_I)
% PLOT_WILSON_COWAN_TRAJECTORY plots stochastic trajectory and WC phase plane.
%
% X       = data matrix (time × variables)
% T       = time vector
% tbin    = bin width for spike counting
% alpha   = leak rate
% w_E     = excitatory weight
% w_I     = inhibitory weight
% h       = external input
% act     = activation function handle, e.g. @(x) 1./(1+exp(-x))

% bin and sum spikes
active = X(:, 2:2:2*N_E);
D = diff(active,1,1);        % diff over T (rows)
D(D < 0) = 0;           % remove negatives
spike_E = [zeros(1,N_E); D];        % prepend first row
[Tb, nspk_E] = timebin_sum(T, spike_E, tbin);
nspk_E = sum(nspk_E, 2);

active = X(:, 2*N_E+2:2:end);
D = diff(active,1,1);        % diff over T (rows)
D(D < 0) = 0;           % remove negatives
spike_I = [zeros(1,N_I); D];        % prepend first row
[Tb, nspk_I] = timebin_sum(T, spike_I, tbin);
nspk_I = sum(nspk_I, 2);


% filtered E/I
E = zeros(size(nspk_E)); 
I = zeros(size(nspk_I));
E(1) = 0.1; 
I(1) = 0;

for t = 2:numel(nspk_E)
    E(t) = (1-alpha*tbin)*E(t-1) + nspk_E(t)/N_E;
    I(t) = (1-alpha*tbin)*I(t-1) + nspk_I(t)/N_I;
end

% plot
f = figure; hold on
plot(E, I, 'g', 'LineWidth',1)

% vector field
[Egrid,Igrid] = meshgrid(linspace(0,1,25));
F = act(w_E*Egrid - w_I*Igrid + h);
dE = -alpha*Egrid + (1 - Egrid).*F;
dI = -alpha*Igrid + (1 - Igrid).*F;
quiver(Egrid,Igrid,dE,dI, 'Color',[0.7 0.7 0.7])

% nullclines
ns = 400;
Evals = linspace(0,1,ns);
Ivals = linspace(0,1,ns);

I_for_Enull = zeros(1,ns);
for k = 1:ns
    fE = act(w_E*Evals(k) - w_I*Ivals + h);
    [~,mi] = min(abs(Evals(k) - fE./(alpha+fE)));
    I_for_Enull(k) = Ivals(mi);
end
plot(Evals, I_for_Enull, 'r', 'LineWidth',2)

E_for_Inull = zeros(1,ns);
for k = 1:ns
    fI = act(w_E*Evals - w_I*Ivals(k) + h);
    [~,mi] = min(abs(Ivals(k) - fI./(alpha+fI)));
    E_for_Inull(k) = Evals(mi);
end
plot(E_for_Inull, Ivals, 'b', 'LineWidth',2)

% deterministic WC
wc = @(t,y)[-alpha*y(1)+(1-y(1))*act(w_E*y(1)-w_I*y(2)+h);
            -alpha*y(2)+(1-y(2))*act(w_E*y(1)-w_I*y(2)+h)];
[tode,yode] = ode45(wc, [0 T(end)], [E(1) I(1)]);
plot(yode(:,1), yode(:,2), 'k--', 'LineWidth',1.5)

xlim([0 1]); ylim([0 1]); axis square
hold off
end
