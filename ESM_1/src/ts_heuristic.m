function H = ts_heuristic(S)
%TS_HEURISTIC  Three-phase heuristic for the turn-optimal Eulerian circuit.
%   1) vertex-wise cheapest pairings (generally several closed trails);
%   2) Kruskal-type merging: repeatedly apply the cheapest re-pairing of two
%      transitions that belong to different closed trails (merge lemma);
%   3) local search over all alternative pairings of one vertex at a time,
%      accepting a change if the result is still one circuit and cheaper.
partner = ts_local_partner(S, ones(S.n,1));
[comp, nc] = ts_components(S, partner);
H.nc0 = nc;
H.delta = zeros(0,1);
while nc > 1
    best = inf; bv = []; bQ = [];
    for v = 1:S.n
        hs = S.inc{v}(:);
        if numel(hs) < 4, continue; end
        T = hs(hs < partner(hs));
        for i = 1:numel(T)-1
            a = T(i); b = partner(a);
            for j = i+1:numel(T)
                c = T(j); d = partner(c);
                if comp(ceil(a/2)) == comp(ceil(c/2)), continue; end
                old = S.Cmin(a,b) + S.Cmin(c,d);
                o1 = S.Cmin(a,c) + S.Cmin(b,d);
                o2 = S.Cmin(a,d) + S.Cmin(b,c);
                if o1 <= o2, dl = o1 - old; Q = [a c; b d];
                else,        dl = o2 - old; Q = [a d; b c]; end
                if dl < best - 1e-9, best = dl; bv = v; bQ = Q; end
            end
        end
    end
    partner = ts_apply(partner, bQ);
    H.delta(end+1,1) = best;
    [comp, nc] = ts_components(S, partner);
end
[cost, ~, ok] = ts_eval(S, partner);
assert(ok, 'merge phase must return a single circuit');
H.costMerge = cost;
H.relaxMerge = ts_mincost(S, partner);
vary = find(S.np > 1)';
improved = true; H.passes = 0; H.moves = 0;
while improved
    improved = false; H.passes = H.passes + 1;
    for v = vary
        for p = 1:S.np(v)
            trial = ts_apply(partner, S.pairs{v}{p});
            [L, F, R] = ts_trace(S, trial, 1);
            if L == S.m && min(F,R) < cost - 1e-9
                partner = trial; cost = min(F,R);
                improved = true; H.moves = H.moves + 1;
            end
        end
    end
end
H.partner = partner;
H.cost = cost;
end
