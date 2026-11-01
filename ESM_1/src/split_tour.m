function R = split_tour(S, G, seq, W, Dw, P, Q)
%SPLIT_TOUR  Optimal capacity-feasible split of an oriented Eulerian circuit
%   (route-first / cluster-second; Beasley 1983, Ulusoy 1985). Each trip
%   leaves the transfer station S, deadheads to the start of its first
%   service street, follows the circuit (services, in-tour deadheads, turns)
%   to the end of its last service street, returns to S and unloads.
%   Every rotation of the circuit is tried; for each, an O(M^2) shortest-path
%   DP over the M service traversals is exact. Times in seconds.
m = numel(seq);
e = ceil(seq/2);
fr = S.hv(seq); to = S.hv(S.oth(seq));
len = G.len(e); svc = G.svc(e);
q = zeros(m,1); q(svc) = W.demand(G.orig(e(svc)));
ttr = len ./ (P.vS*svc + P.vD*(~svc));
trn = zeros(m,1);
for k = 1:m, trn(k) = S.C(S.oth(seq(k)), seq(mod(k,m)+1)); end
% station <-> vertex deadhead (time and distance)
ng = numel(W.gate);
tS = inf(G.n,1); dS = inf(G.n,1);
for g = 1:ng
    t = Dw(:, W.gate(g))/P.vD + W.gateLen(g)/P.vA;
    better = t < tS;
    tS(better) = t(better);
    dS(better) = Dw(better, W.gate(g)) + W.gateLen(g);
end
% doubled arrays for cyclic rotation
ttr2 = [ttr; ttr]; trn2 = [trn; trn]; q2 = [q; q]; svc2 = [svc; svc];
len2 = [len; len]; fr2 = [fr; fr]; to2 = [to; to];
A = [0; cumsum(ttr2)]; Bt = [0; cumsum(trn2)];
Ld = [0; cumsum(len2 .* ~svc2)]; Ls = [0; cumsum(len2 .* svc2)];
sp = find(svc);
best = inf;
for r = sp'
    ps = r - 1 + find(svc2(r:r+m-1));        % service positions in this rotation
    M = numel(ps);
    Qc = [0; cumsum(q2(ps))];
    f = [0; inf(M,1)]; arg = zeros(M,1);
    for j = 1:M
        for i = j:-1:1
            load = Qc(j+1) - Qc(i);
            if load > Q + 1e-9, break; end
            t = tS(fr2(ps(i))) + (A(ps(j)+1) - A(ps(i))) + (Bt(ps(j)) - Bt(ps(i))) ...
                + tS(to2(ps(j))) + P.tUnload;
            if f(i) + t < f(j+1) - 1e-9, f(j+1) = f(i) + t; arg(j) = i; end
        end
    end
    if f(M+1) < best - 1e-9
        best = f(M+1); bestR = r; bestPs = ps; bestArg = arg; bestQc = Qc;
    end
end
% reconstruct the trips of the best rotation
ps = bestPs; M = numel(ps); j = M; trips = zeros(0,2);
while j > 0, i = bestArg(j); trips = [i j; trips]; j = i - 1; end %#ok<AGROW>
nt = size(trips,1);
R.rotation = bestR; R.total = best; R.ntrips = nt;
R.trip = struct('first', {}, 'last', {}, 'load', {}, 'time', {}, 'pos', {});
R.tService = 0; R.tDeadIn = 0; R.tAccess = 0; R.tTurn = 0; R.tUnload = nt*P.tUnload;
R.dService = 0; R.dDeadIn = 0; R.dAccess = 0;
for t = 1:nt
    i = trips(t,1); j = trips(t,2);
    a = ps(i); b = ps(j);
    seg = a:b;
    R.trip(t).first = i; R.trip(t).last = j;
    R.trip(t).load = bestQc(j+1) - bestQc(i);
    R.trip(t).pos = mod(seg-1, m) + 1;     % circuit positions covered
    tin = A(b+1) - A(a); tt = Bt(b) - Bt(a);
    tacc = tS(fr2(a)) + tS(to2(b));
    R.trip(t).time = tin + tt + tacc + P.tUnload;
    R.trip(t).access = [dS(fr2(a)) dS(to2(b))];
    ts_ = sum(ttr2(seg) .* svc2(seg));
    R.tService = R.tService + ts_;
    R.tDeadIn = R.tDeadIn + (tin - ts_);
    R.tTurn = R.tTurn + tt; R.tAccess = R.tAccess + tacc;
    R.dService = R.dService + (Ls(b+1) - Ls(a));
    R.dDeadIn = R.dDeadIn + (Ld(b+1) - Ld(a));
    R.dAccess = R.dAccess + dS(fr2(a)) + dS(to2(b));
end
R.dTotal = R.dService + R.dDeadIn + R.dAccess;
end
