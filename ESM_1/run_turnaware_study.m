function R = run_turnaware_study(doBench)
%RUN_TURNAWARE_STUDY  Reproduces every number, table and figure of the article
%   "Selecting among equally short collection routes: turn-penalty-optimal
%    Eulerian circuits via transition systems, with capacity-feasible tour
%    splitting, for kerbside waste collection".
%
%   R = run_turnaware_study        full run (about 3-6 min in GNU Octave)
%   R = run_turnaware_study(false) skips the random-instance benchmark
%
%   Runs unchanged in MATLAB R2019b+ (exact solver: intlinprog) and in
%   GNU Octave 8+ (exact solver: glpk). All randomness comes from the
%   Wichmann-Hill generator in src/wh_rand.m, so every number except
%   wall-clock time is identical on both platforms.
%
%   Outputs: results/run_log.txt, results/*.csv, results/results.mat,
%            figs/fig01..fig12 (.pdf and .png).

if nargin < 1, doBench = true; end
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, 'src'));
rdir = fullfile(here, 'results'); fdir = fullfile(here, 'figs');
if ~exist(rdir, 'dir'), mkdir(rdir); end
if ~exist(fdir, 'dir'), mkdir(fdir); end
logmsg('#open', fullfile(rdir, 'run_log.txt'));
isOct = exist('OCTAVE_VERSION', 'builtin') > 0;
if isOct, env = ['GNU Octave ' OCTAVE_VERSION]; else, env = ['MATLAB ' version]; end
logmsg('Turn-penalty-optimal Eulerian collection routes\n');
logmsg('================================================\n');
logmsg('Environment : %s\nRun started : %s\n\n', env, datestr(now, 'yyyy-mm-dd HH:MM:SS'));

%% ------------------------------------------------------------------ params
P = struct('tS', 0, 'tL', 6, 'tR', 18, 'tU', 40, 'degS', 30, 'degU', 150);
P.vS = 6/3.6;     % collection (service) speed, m/s
P.vD = 20/3.6;    % in-ward deadhead speed, m/s
P.vA = 30/3.6;    % access-road speed to the transfer station, m/s
P.tUnload = 900;  % unloading time per trip, s
Qcap = 1500;      % truck payload, kg
Nsample = 5000;   % uniform random Eulerian circuits
Nsplit = 400;     % of which split into trips for the time distribution

%% ------------------------------------------------- 1. network and parity
W = ward_network(); n = size(W.xy,1); m = size(W.E,1);
deg = accumarray(W.E(:), 1, [n 1]);
odd = find(mod(deg,2))';
logmsg('1. Ward network\n---------------\n');
logmsg('n = %d intersections, m = %d streets, total length %d m, demand %.1f kg\n', ...
    n, m, sum(W.len), sum(W.demand));
logmsg('  ID   Street    Length  Class                  Demand(kg)\n');
for e = 1:m
    logmsg('  e%-3d %-4s-%-4s %5d   %-22s %7.2f\n', e, W.name{W.E(e,1)}, ...
        W.name{W.E(e,2)}, W.len(e), W.className{W.class(e)}, W.demand(e));
end
dg = [W.name; num2cell(deg')];
logmsg('Degrees: %s\n', sprintf('%s:%d ', dg{:}));
logmsg('Handshake: sum deg = %d = 2m = %d -> %s\n', sum(deg), 2*m, pf(sum(deg)==2*m));
logmsg('Odd vertices (%d): %s\n\n', numel(odd), strjoin(W.name(odd), ', '));
for c = 1:3
    logmsg('  class %d (%s): %d streets, %d m, %.1f kg\n', c, W.className{c}, ...
        sum(W.class==c), sum(W.len(W.class==c)), sum(W.demand(W.class==c)));
end
R.W = W; R.deg = deg; R.odd = odd; R.P = P;

%% ------------------------------------------------- 2. parity repair (T-join)
[D, Nx] = apsp_floyd(n, W.E, W.len);
Mall = all_matchings(odd);
nm = numel(Mall);
mw = zeros(nm,1);
for k = 1:nm, mw(k) = sum(D(sub2ind([n n], Mall{k}(:,1), Mall{k}(:,2)))); end
[pr, dh] = min_matching_dp(odd, D);
[~, kmin] = min(mw);
logmsg('\n2. Closed parity repair\n-----------------------\n');
logmsg('%d perfect matchings of the odd set; weights %d .. %d m\n', nm, min(mw), max(mw));
logmsg('Subset-DP optimum %d m, enumeration optimum %d m -> %s\n', dh, mw(kmin), pf(abs(dh-mw(kmin))<1e-9));
for r = 1:size(pr,1)
    p = sp_path(Nx, pr(r,1), pr(r,2));
    logmsg('  %s-%s : %d m via %s\n', W.name{pr(r,1)}, W.name{pr(r,2)}, D(pr(r,1),pr(r,2)), ...
        strjoin(W.name(p), '-'));
end
G = eulerize_closed(n, W.E, W.len, pr, Nx);
deg2 = accumarray(G.E(:), 1, [n 1]);
logmsg('Duplicated streets (%d): %s ; deadhead %d m\n', numel(G.dup), ...
    strjoin(arrayfun(@(e) sprintf('e%d', e), G.dup', 'UniformOutput', false), ' '), sum(G.len(~G.svc)));
logmsg('Repaired multigraph: %d edges, %d m; all degrees even -> %s\n', G.m, sum(G.len), pf(all(mod(deg2,2)==0)));
R.matchW = mw; R.G = G; R.deg2 = deg2;

%% ------------------------------------------------- 3. transition systems
S = ts_prepare(G, W.xy, P);
logmsg('\n3. Turn model and transition systems\n------------------------------------\n');
logmsg('Penalties (s): straight %d, left %d, right %d, U-turn %d; |theta|<=%d straight, >%d U\n', ...
    P.tS, P.tL, P.tR, P.tU, P.degS, P.degU);
logmsg('Pairings per vertex: %s\n', mat2str(S.np'));
logmsg('Transition systems: prod (d-1)!! = %.4g (log10 = %.3f)\n', 10^S.logNTS, S.logNTS);
logmsg('Local lower bound LB = %d s\n', S.LB);
R.S = S;

%% ------------------------------------------------- 4. baseline / heuristic / exact
pb = hierholzer_partner(S, G, 1);
[cBase, seqBase, okB] = ts_eval(S, pb);
Hh = ts_heuristic(S);
X = ts_ilp(S, G);
[cOpt, seqOpt, okO] = ts_eval(S, X.partner);
logmsg('\n4. Circuits\n-----------\n');
logmsg('Hierholzer baseline : %d s (single circuit %s)\n', cBase, pf(okB));
logmsg('Heuristic           : start %d closed trails, %d merges, after merge %d s, after LS %d s (%d moves, %d passes)\n', ...
    Hh.nc0, numel(Hh.delta), Hh.costMerge, Hh.cost, Hh.moves, Hh.passes);
logmsg('Exact 0-1 programme : %d s, %d variables (%d transitions), %d rounds, %d cuts, %.2f s -> ok %s\n', ...
    X.cost, X.nvar, X.ntrans, X.rounds, X.cuts, X.time, pf(X.ok));
logmsg('Re-evaluated optimum: %d s (single circuit %s, equals ILP objective %s)\n', cOpt, pf(okO), pf(abs(cOpt-X.cost)<1e-9));
logmsg('Gap to LB: optimum %.1f%% above LB; baseline %.1f%% above optimum\n', ...
    100*(cOpt/S.LB-1), 100*(cBase/cOpt-1));
[clsB, cstB, angB] = circuit_turns(S, seqBase);
[clsO, cstO, angO] = circuit_turns(S, seqOpt);
[~, seqH] = ts_eval(S, Hh.partner); clsH = circuit_turns(S, seqH);
cnt = @(c) [sum(c==1) sum(c==2) sum(c==3) sum(c==4)];
logmsg('Turn counts [S L R U]: baseline %s, heuristic %s, optimum %s\n', ...
    mat2str(cnt(clsB)), mat2str(cnt(clsH)), mat2str(cnt(clsO)));
% validation of the optimum
V = validate_circuit(S, G, seqOpt, W);
logmsg('Validation of optimum:\n');
for k = 1:size(V,1), logmsg('  %-52s %-10s %-10s %s\n', V{k,:}); end
R.base = struct('cost', cBase, 'seq', seqBase, 'cls', clsB, 'cst', cstB, 'ang', angB);
R.opt  = struct('cost', cOpt, 'seq', seqOpt, 'cls', clsO, 'cst', cstO, 'ang', angO);
R.heur = Hh; R.ilp = rmfield(X, {'x'}); R.validation = V;
logmsg('Optimal circuit (from %s): %s\n', W.name{S.hv(seqOpt(1))}, ...
    strjoin(W.name(S.hv([seqOpt; seqOpt(1)])), '-'));
logmsg('Turn sequence: %s\n', turnstr(clsO));

%% ------------------------------------------------- 5. uniform sampling
tic5 = tic;
Ms = ts_sample(S, Nsample, 50*Nsample, 2026);
nEC = Ms.p * 10^S.logNTS; nECci = Ms.pCI * 10^S.logNTS;
logmsg('\n5. Uniform random Eulerian circuits\n-----------------------------------\n');
logmsg('%d accepted of %d draws (p = %.4f, 95%% CI %.4f-%.4f), %.1f s\n', Ms.accepted, Ms.trials, ...
    Ms.p, Ms.pCI(1), Ms.pCI(2), toc(tic5));
logmsg('Estimated number of Eulerian circuits ec(G) = %.3g (95%% CI %.3g-%.3g)\n', nEC, nECci(1), nECci(2));
logmsg('Turn delay over random circuits: mean %.1f s, sd %.1f, median %g, min %g, max %g\n', ...
    mean(Ms.cost), std(Ms.cost), median(Ms.cost), min(Ms.cost), max(Ms.cost));
logmsg('Fraction of random circuits within 10%% of optimum: %.4f ; optimum hit: %d of %d\n', ...
    mean(Ms.cost <= 1.1*cOpt), sum(Ms.cost == cOpt), Ms.accepted);
logmsg('Optimum percentile among random circuits: %.3f%%; baseline percentile: %.1f%%\n', ...
    100*mean(Ms.cost < cOpt), 100*mean(Ms.cost < cBase));
logmsg('Mean turn counts [S L R U] of random circuits: %s\n', mat2str(round(10*[mean(Ms.nS) mean(Ms.nL) mean(Ms.nR) mean(Ms.nU)])/10));
R.sample = rmfield(Ms, 'partner'); R.nEC = nEC; R.nECci = nECci;

%% ------------------------------------------------- 6. joint repair analysis
logmsg('\n6. All %d repairs: deadhead versus optimal turn delay\n------------------------------------------------\n', nm);
JA = zeros(nm, 6);
for k = 1:nm
    Gk = eulerize_closed(n, W.E, W.len, Mall{k}, Nx);
    Sk = ts_prepare(Gk, W.xy, P);
    Xk = ts_ilp(Sk, Gk);
    pbk = hierholzer_partner(Sk, Gk, 1); cbk = ts_eval(Sk, pbk);
    dk = sum(Gk.len(~Gk.svc));
    tClosed = sum(W.len)/P.vS + dk/P.vD + Xk.cost;     % closed-tour time, s
    JA(k,:) = [k, dk, Xk.cost, cbk, tClosed, Gk.m];
    logmsg('  M%-2d %-28s deadhead %4d m  turn* %3d s  (Hierholzer %3d s)  tour %.1f min\n', k, ...
        strjoin(arrayfun(@(r) [W.name{Mall{k}(r,1)} W.name{Mall{k}(r,2)}], 1:size(Mall{k},1), 'UniformOutput', false), ','), ...
        dk, Xk.cost, cbk, tClosed/60);
end
[~, kT] = min(JA(:,5));
logmsg('Distance-optimal repair M%d: tour %.2f min; time-optimal repair M%d: tour %.2f min\n', ...
    kmin, JA(kmin,5)/60, kT, JA(kT,5)/60);
R.joint = JA; R.Mall = Mall; R.kmin = kmin; R.kT = kT;

%% ------------------------------------------------- 7. capacity split
Sb = split_tour(S, G, seqBase, W, D, P, Qcap);
So = split_tour(S, G, seqOpt, W, D, P, Qcap);
logmsg('\n7. Capacity-feasible trips (Q = %d kg, total demand %.1f kg)\n---------------------------------------------\n', Qcap, sum(W.demand));
report_split('Baseline', Sb);
report_split('Optimum ', So);
tic7 = tic;
rs = zeros(Nsplit, 7);
for k = 1:Nsplit
    [~, sk] = ts_eval(S, Ms.partner(:,k));
    Sk = split_tour(S, G, sk, W, D, P, Qcap);
    rs(k,:) = [Sk.total, Sk.tService, Sk.tDeadIn, Sk.tAccess, Sk.tTurn, Sk.tUnload, Sk.ntrips];
end
logmsg('Random circuits (%d): shift %.2f +- %.2f min (min %.2f, max %.2f), turn time in trips %.1f s; %.1f s\n', ...
    Nsplit, mean(rs(:,1))/60, std(rs(:,1))/60, min(rs(:,1))/60, max(rs(:,1))/60, mean(rs(:,5)), toc(tic7));
logmsg('Shift saving of optimum: %.1f s vs baseline (%.2f%%), %.1f s vs random mean (%.2f%%)\n', ...
    Sb.total - So.total, 100*(1 - So.total/Sb.total), mean(rs(:,1)) - So.total, 100*(1 - So.total/mean(rs(:,1))));
logmsg('Annual (300 collection days): %.1f h saved vs baseline, %.1f h vs random mean\n', ...
    300*(Sb.total - So.total)/3600, 300*(mean(rs(:,1)) - So.total)/3600);
logmsg('Right turns avoided per year vs baseline: %d\n', 300*(sum(clsB==3) - sum(clsO==3)));
R.splitBase = Sb; R.splitOpt = So; R.splitRand = rs; R.Q = Qcap;

%% ------------------------------------------------- 8. penalty sensitivity
tLs = [0 3 6 9 12]; tRs = [6 12 18 24 30 36];
SR = zeros(numel(tLs), numel(tRs)); NR = SR; NRb = SR;
for i = 1:numel(tLs)
    for j = 1:numel(tRs)
        Pk = P; Pk.tL = tLs(i); Pk.tR = tRs(j);
        Sk = ts_prepare(G, W.xy, Pk);
        Xk = ts_ilp(Sk, G);
        [~, sk] = ts_eval(Sk, Xk.partner); ck = circuit_turns(Sk, sk);
        [cb_, sbk] = ts_eval(Sk, pb); cbk = circuit_turns(Sk, sbk);
        SR(i,j) = 100*(1 - Xk.cost/cb_); NR(i,j) = sum(ck==3); NRb(i,j) = sum(cbk==3);
    end
end
logmsg('\n8. Sensitivity to turn penalties (reduction %% vs baseline / right turns in optimum)\n');
logmsg('   tL\\tR %s\n', sprintf('%7d', tRs));
for i = 1:numel(tLs)
    logmsg('   %4d  %s\n', tLs(i), sprintf('%5.1f/%d', [SR(i,:); NR(i,:)]));
end
R.sens = struct('tL', tLs, 'tR', tRs, 'red', SR, 'nR', NR, 'nRbase', NRb);

%% ------------------------------------------------- 9. capacity sensitivity
Qs = 800:100:3600;
QS = zeros(numel(Qs), 4);
for k = 1:numel(Qs)
    a = split_tour(S, G, seqBase, W, D, P, Qs(k));
    b = split_tour(S, G, seqOpt, W, D, P, Qs(k));
    QS(k,:) = [a.total b.total a.ntrips b.ntrips];
end
logmsg('\n9. Capacity sensitivity: Q = %d..%d kg; trips %d..%d; saving %.1f..%.1f s\n', ...
    Qs(1), Qs(end), min(QS(:,4)), max(QS(:,4)), min(QS(:,1)-QS(:,2)), max(QS(:,1)-QS(:,2)));
% split-aware local search wherever the turn-optimal circuit loses to the
% baseline, and at the design payload
QA = QS(:,2); TA = QS(:,4); NA = zeros(numel(Qs),1);
for k = find(QS(:,2) > QS(:,1) | Qs(:) == Qcap)'
    [~, Sk, nmv] = split_ls(S, G, X.partner, W, D, P, Qs(k));
    QA(k) = min(QS(k,2), Sk.total); TA(k) = Sk.ntrips; NA(k) = nmv;
    logmsg('  Q = %d kg: baseline %.2f min (%d trips), turn-optimal %.2f min (%d trips), split-aware LS %.2f min (%d trips, %d moves)\n', ...
        Qs(k), QS(k,1)/60, QS(k,3), QS(k,2)/60, QS(k,4), Sk.total/60, Sk.ntrips, nmv);
    [clsA, cstA] = circuit_turns(S, Sk.seq);
    logmsg('      split-aware circuit: turn delay %d s, turns [S L R U] = %s\n', sum(cstA), mat2str(histc(clsA', 1:4)));
    if Qs(k) == Qcap, R.splitAware = Sk; R.splitAwareTurn = [sum(cstA) histc(clsA', 1:4)]; end
end
logmsg('After split-aware LS: saving vs baseline %.1f..%.1f s over all payloads\n', min(QS(:,1)-QA), max(QS(:,1)-QA));
R.capsens = [Qs(:) QS QA TA NA];

%% ------------------------------------------------- 10. benchmark
if doBench
    sz = [3 4; 4 5; 5 6; 6 7]; nrep = 10;
    BM = zeros(0, 14);
    logmsg('\n10. Random-instance benchmark (%d instances per size, diagonal prob 0.35)\n', nrep);
    for s = 1:size(sz,1)
        for rep = 1:nrep
            Wr = random_ward(sz(s,1), sz(s,2), 0.35, 1000*s + rep);
            o = solve_instance(Wr, P);
            bb = NaN;
            if o.logNTS <= 6.5
                B = ts_bnb(o.S, o.Xpartner, inf, 2e6);
                if B.proven, bb = B.cost; end
            end
            BM(end+1,:) = [s rep o.n o.m o.m2 o.nodd o.logNTS o.LB o.base o.heur o.exact o.tILP o.tHeur bb]; %#ok<AGROW>
        end
        b = BM(BM(:,1)==s,:);
        logmsg('  %dx%d: n=%d, m=%.1f, logTS=%.1f | gap heur %.1f%% (max %.1f), baseline +%.1f%%, LB gap %.1f%% | ILP %.2f s (max %.2f), heur %.2f s | B&B checks %d/%d agree\n', ...
            sz(s,1), sz(s,2), b(1,3), mean(b(:,4)), mean(b(:,7)), 100*mean(b(:,10)./b(:,11)-1), ...
            100*max(b(:,10)./b(:,11)-1), 100*mean(b(:,9)./b(:,11)-1), 100*mean(1-b(:,8)./b(:,11)), ...
            mean(b(:,12)), max(b(:,12)), mean(b(:,13)), sum(b(:,14)==b(:,11)), sum(~isnan(b(:,14))));
    end
    R.bench = BM; R.benchSizes = sz;
end

%% ------------------------------------------------- tables
write_csv(fullfile(rdir, 'table_streets.csv'), {'id','from','to','length_m','class','demand_kg','duplicated'}, ...
    [(1:m)' W.E W.len W.class W.demand ismember((1:m)', G.dup)]);
write_csv(fullfile(rdir, 'table_repairs.csv'), {'matching','deadhead_m','turn_opt_s','turn_hierholzer_s','closed_tour_s','edges'}, JA);
write_csv(fullfile(rdir, 'table_optimal_circuit.csv'), {'step','edge','from','to','length_m','service','turn_class','turn_s','angle_deg'}, ...
    [(1:G.m)' ceil(seqOpt/2) S.hv(seqOpt) S.hv(S.oth(seqOpt)) G.len(ceil(seqOpt/2)) G.svc(ceil(seqOpt/2)) clsO cstO round(angO)]);
write_csv(fullfile(rdir, 'table_sample.csv'), {'cost_s','nS','nL','nR','nU'}, [Ms.cost Ms.nS Ms.nL Ms.nR Ms.nU]);
write_csv(fullfile(rdir, 'table_capacity.csv'), {'Q_kg','shift_base_s','shift_opt_s','trips_base','trips_opt','shift_splitaware_s','trips_splitaware','ls_moves'}, R.capsens);
tt = zeros(0,6);
for t = 1:So.ntrips, tt(end+1,:) = [2 t So.trip(t).load So.trip(t).time numel(So.trip(t).pos) So.trip(t).first]; end %#ok<AGROW>
for t = 1:Sb.ntrips, tt(end+1,:) = [1 t Sb.trip(t).load Sb.trip(t).time numel(Sb.trip(t).pos) Sb.trip(t).first]; end %#ok<AGROW>
write_csv(fullfile(rdir, 'table_trips.csv'), {'circuit(1=base,2=opt)','trip','load_kg','time_s','traversals','first_service'}, tt);
if doBench
    write_csv(fullfile(rdir, 'table_benchmark.csv'), {'size','rep','n','m','m2','nodd','log10TS','LB','base','heur','exact','tILP','tHeur','bnb'}, BM);
end

%% ------------------------------------------------- figures
make_figures(R, W, G, S, P, Ms, Sb, So, fdir, doBench);

%% ------------------------------------------------- summary
fs = fopen(fullfile(rdir, 'summary.txt'), 'w');
fprintf(fs, 'Environment              : %s\n', env);
fprintf(fs, 'Network                  : n=%d, m=%d, %d m, %d odd vertices\n', n, m, sum(W.len), numel(odd));
fprintf(fs, 'Closed repair            : +%d m deadhead (%d duplicated streets)\n', dh, numel(G.dup));
fprintf(fs, 'Transition systems       : 10^%.3f ; ec(G) ~ %.3g\n', S.logNTS, nEC);
fprintf(fs, 'Turn delay [s]           : LB %d | optimum %d | heuristic %d | Hierholzer %d | random mean %.1f\n', S.LB, cOpt, Hh.cost, cBase, mean(Ms.cost));
fprintf(fs, 'Right turns              : optimum %d | Hierholzer %d | random mean %.1f\n', sum(clsO==3), sum(clsB==3), mean(Ms.nR));
fprintf(fs, 'Shift time [min]         : optimum %.2f | Hierholzer %.2f | random mean %.2f (%d trips)\n', So.total/60, Sb.total/60, mean(rs(:,1))/60, So.ntrips);
fclose(fs);
save(fullfile(rdir, 'results.mat'), 'R');
logmsg('\nDone: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
logmsg('#close');
end

%% =================================================================== helpers
function s = pf(tf)
if tf, s = 'PASS'; else, s = 'FAIL'; end
end

function report_split(lbl, Sx)
logmsg('%s: %d trips, shift %.2f min = service %.2f + in-ward deadhead %.2f + access %.2f + turns %.2f + unloading %.2f\n', ...
    lbl, Sx.ntrips, Sx.total/60, Sx.tService/60, Sx.tDeadIn/60, Sx.tAccess/60, Sx.tTurn/60, Sx.tUnload/60);
logmsg('          distance %.0f m = service %.0f + in-ward deadhead %.0f + access %.0f; loads %s kg\n', ...
    Sx.dTotal, Sx.dService, Sx.dDeadIn, Sx.dAccess, mat2str(round([Sx.trip.load])));
end

function V = validate_circuit(S, G, seq, W)
e = ceil(seq/2);
cntE = accumarray(e(:), 1, [G.m 1]);
contig = all(S.hv(S.oth(seq)) == S.hv(seq([2:end 1])));
streetCnt = accumarray(G.orig(e), 1, [size(W.E,1) 1]);
need = 1 + accumarray(G.dup(:), 1, [size(W.E,1) 1]);
svcCnt = accumarray(G.orig(e(G.svc(e))), 1, [size(W.E,1) 1]);
V = {
 'Traversals = edges of repaired multigraph', sprintf('%d', G.m), sprintf('%d', numel(seq)), pf(numel(seq)==G.m);
 'Every multigraph edge used exactly once', 'yes', yn(all(cntE==1)), pf(all(cntE==1));
 'Consecutive traversals share a vertex (closed walk)', 'yes', yn(contig), pf(contig);
 'Each street driven 1 + (#duplicates) times', 'yes', yn(all(streetCnt==need)), pf(all(streetCnt==need));
 'Each street serviced exactly once', 'yes', yn(all(svcCnt==1)), pf(all(svcCnt==1));
 'Route length (m) = streets + T-join', sprintf('%d', sum(W.len)+sum(G.len(~G.svc))), sprintf('%d', sum(G.len(e))), ...
     pf(sum(G.len(e)) == sum(W.len)+sum(G.len(~G.svc)))};
end

function s = turnstr(c)
L = 'SLRU';
s = L(c(:)');
end

function s = yn(tf)
if tf, s = 'yes'; else, s = 'no'; end
end

function write_csv(fn, hdr, M)
fid = fopen(fn, 'w');
fprintf(fid, '%s\n', strjoin(hdr, ','));
for i = 1:size(M,1)
    fprintf(fid, '%s\n', strjoin(arrayfun(@(x) num2str(x, 10), M(i,:), 'UniformOutput', false), ','));
end
fclose(fid);
end
