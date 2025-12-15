function [T, X] = gillespie(reactions, stoich, x0, t_end)

rng('shuffle');

t = 0;
x = x0(:);

Tcell = {t};
Xcell = {x.'};

k = 1;
while t < t_end
    % propensities
    a = zeros(size(stoich,2),1);
    for j = 1:numel(a)
        a(j) = reactions{j}(x);
    end
    a0 = sum(a);
    if a0 == 0, break; end

    % Gillespie step
    tau = -log(rand()) / a0;
    k_rxn = find(rand()*a0 <= cumsum(a), 1);

    x = x + stoich(:, k_rxn);
    t = t + tau;

    % store
    k = k + 1;
    Tcell{k,1} = t;
    Xcell{k,1} = x.';

    msg = sprintf('\rTime Progress (ms): %f / %f', t, t_end); 
    fprintf([repmat('\b',1,100) msg]);
end
fprintf('\n');

% ---- concatenate once ----
T = vertcat(Tcell{:});   % (N,1)
X = vertcat(Xcell{:});   % (N,numSpecies)

end