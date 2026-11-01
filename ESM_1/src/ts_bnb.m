function B = ts_bnb(S, partner0, ub0, maxNodes)
%TS_BNB  Exact branch-and-bound over transition systems for the minimum
%   turn-cost Eulerian circuit. Vertices are branched in decreasing order of
%   their number of pairings; pairings are tried in increasing order of their
%   relaxed cost, so that one failed bound test closes the whole level.
%   Bound: relaxed cost of fixed vertices + sum of vertex-wise minima of the
%   free vertices (Proposition 4.3). Leaves are evaluated exactly.
if nargin < 4, maxNodes = inf; end
tic0 = tic;
vary = find(S.np > 1);
[~, o] = sort(-S.np(vary)); order = vary(o);
fixed = find(S.np == 1);
K = numel(order);
lbv = zeros(K,1);
for k = 1:K, lbv(k) = S.mc{order(k)}(1); end
suffix = [flipud(cumsum(flipud(lbv))); 0];
base = 0;
for v = fixed', base = base + S.mc{v}(1); end
partner = partner0;
for v = fixed', partner = ts_apply(partner, S.pairs{v}{1}); end
best = ub0; bestP = partner0;
tol = 1e-9;
ch = zeros(K,1); partial = zeros(K+1,1); partial(1) = base;
nodes = 0; leaves = 0; circuits = 0;
hist = [0 best];
proven = true;
k = 1;
while k >= 1
    ch(k) = ch(k) + 1;
    v = order(k);
    if ch(k) > S.np(v)
        ch(k) = 0; k = k - 1; continue;
    end
    c = partial(k) + S.mc{v}(ch(k));
    if c + suffix(k+1) >= best - tol
        ch(k) = 0; k = k - 1; continue;      % sorted order: prune the level
    end
    partner = ts_apply(partner, S.pairs{v}{ch(k)});
    nodes = nodes + 1;
    if nodes > maxNodes, proven = false; break; end
    if k == K
        leaves = leaves + 1;
        [L, F, R] = ts_trace(S, partner, 1);
        if L == S.m
            circuits = circuits + 1;
            val = min(F, R);
            if val < best - tol
                best = val; bestP = partner;
                hist(end+1,:) = [nodes best]; %#ok<AGROW>
            end
        end
    else
        partial(k+1) = c;
        k = k + 1;
    end
end
B.cost = best; B.partner = bestP; B.nodes = nodes; B.leaves = leaves;
B.circuitLeaves = circuits; B.hist = hist; B.proven = proven;
B.time = toc(tic0); B.K = K;
end
