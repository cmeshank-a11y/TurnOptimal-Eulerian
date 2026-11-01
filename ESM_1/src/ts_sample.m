function M = ts_sample(S, N, maxTrials, seed)
%TS_SAMPLE  Uniform random Eulerian circuits by rejection: draw every vertex
%   pairing independently and uniformly, keep the draw iff it is a single
%   closed trail. Accepted draws are uniform over Eulerian circuits
%   (Proposition 5.2); the acceptance rate estimates ec(G)/prod((d-1)!!).
wh_rand('seed', seed);
vary = find(S.np > 1);
M.cost = zeros(N,1); M.nR = zeros(N,1); M.nL = zeros(N,1);
M.nS = zeros(N,1); M.nU = zeros(N,1); M.partner = zeros(S.H, N);
acc = 0; trials = 0;
while acc < N && trials < maxTrials
    trials = trials + 1;
    idx = ones(S.n,1);
    idx(vary) = floor(wh_rand(numel(vary)) .* S.np(vary)) + 1;
    partner = ts_local_partner(S, idx);
    [L, F, R] = ts_trace(S, partner, 1);
    if L == S.m
        acc = acc + 1;
        [c, seq] = ts_eval(S, partner);
        cls = circuit_turns(S, seq);
        M.cost(acc) = c;
        M.nS(acc) = sum(cls==1); M.nL(acc) = sum(cls==2);
        M.nR(acc) = sum(cls==3); M.nU(acc) = sum(cls==4);
        M.partner(:,acc) = partner;
    end
end
f = {'cost','nR','nL','nS','nU'};
for i = 1:numel(f), M.(f{i}) = M.(f{i})(1:acc); end
M.partner = M.partner(:,1:acc);
M.accepted = acc; M.trials = trials;
M.p = acc / trials;
M.pCI = M.p + [-1 1] * 1.96 * sqrt(M.p*(1-M.p)/trials);
end
