function out = solve_instance(W, P)
%SOLVE_INSTANCE  Closed parity repair + baseline / heuristic / exact turn
%   optimisation for one network W. Used by the benchmark and sensitivity runs.
n = size(W.xy,1);
deg = accumarray(W.E(:), 1, [n 1]);
odd = find(mod(deg,2));
[D, Nx] = apsp_floyd(n, W.E, W.len);
[pr, dh] = min_matching_dp(odd(:)', D);
G = eulerize_closed(n, W.E, W.len, pr, Nx);
S = ts_prepare(G, W.xy, P);
pb = hierholzer_partner(S, G, 1);
out.base = ts_eval(S, pb);
t = tic; Hh = ts_heuristic(S); out.tHeur = toc(t);
out.heur = Hh.cost;
X = ts_ilp(S, G);
out.exact = X.cost; out.ok = X.ok; out.tILP = X.time;
out.rounds = X.rounds; out.cuts = X.cuts;
out.LB = S.LB; out.n = n; out.m = size(W.E,1); out.m2 = G.m;
out.nodd = numel(odd); out.deadhead = dh; out.logNTS = S.logNTS;
out.S = S; out.G = G; out.Xpartner = X.partner;
end
