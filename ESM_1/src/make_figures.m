function make_figures(R, W, G, S, P, Ms, Sb, So, fdir, doBench)
%MAKE_FIGURES  Twelve figures of the article (PDF + PNG) from the results R.
C.blue = [42 120 214]/255; C.orange = [235 104 52]/255; C.aqua = [27 175 122]/255;
C.yellow = [237 161 0]/255; C.violet = [98 80 214]/255; C.red = [227 73 72]/255;
C.grey = [0.62 0.62 0.60]; C.dark = [0.20 0.20 0.20]; C.light = [0.86 0.86 0.84];
clsCol = [C.blue; C.orange; C.aqua];          % street classes
turnCol = [0.55 0.55 0.55; C.blue; C.red; C.violet];   % S L R U
turnName = {'Straight', 'Left', 'Right', 'U-turn'};
set(0, 'defaultaxesfontsize', 9, 'defaulttextfontsize', 9);
n = size(W.xy,1); m = size(W.E,1);
xy = W.xy;

%% fig01 ward network
f = newfig(6.4, 5.0); hold on; axis equal off;
for e = 1:m
    a = W.E(e,1); b = W.E(e,2);
    plot(xy([a b],1), xy([a b],2), '-', 'color', clsCol(W.class(e),:), 'linewidth', 2.4);
    mid = mean(xy([a b],:)); nv = perp(xy(a,:), xy(b,:));
    text(mid(1)+11*nv(1), mid(2)+11*nv(2), sprintf('%d', W.len(e)), 'fontsize', 6.5, ...
        'color', [0.35 0.35 0.35], 'horizontalalignment', 'center');
end
for g = 1:numel(W.gate)
    v = W.gate(g);
    plot([xy(v,1) xy(v,1)], [xy(v,2) xy(v,2)-70], '--', 'color', C.dark, 'linewidth', 1);
    text(xy(v,1)+6, xy(v,2)-80, sprintf('to S: %d m', W.gateLen(g)), 'fontsize', 7, 'color', C.dark);
end
drawnodes(xy, W.name, R.odd, C);
h1 = plot(nan, nan, '-', 'color', clsCol(1,:), 'linewidth', 2.4);
h2 = plot(nan, nan, '-', 'color', clsCol(2,:), 'linewidth', 2.4);
h3 = plot(nan, nan, '-', 'color', clsCol(3,:), 'linewidth', 2.4);
h4 = plot(nan, nan, 'o', 'markerfacecolor', C.red, 'markeredgecolor', 'k', 'markersize', 7);
legend([h1 h2 h3 h4], {'Residential (medium)', 'Residential (high)', 'Commercial', 'Odd-degree intersection'}, ...
    'location', 'northoutside', 'orientation', 'horizontal', 'fontsize', 7.5);
legend boxoff
xlim([-60 570]); ylim([-110 500]);
savefig(f, fdir, 'fig01_ward_network');

%% fig02 turn model
f = newfig(7.0, 2.8);
subplot(1,2,1); hold on; box on;
th = -180:0.5:180; pen = zeros(size(th));
for k = 1:numel(th)
    t = th(k);
    if abs(t) <= P.degS, pen(k) = P.tS; elseif abs(t) > P.degU, pen(k) = P.tU;
    elseif t > 0, pen(k) = P.tL; else, pen(k) = P.tR; end
end
plot(th, pen, '-', 'color', C.dark, 'linewidth', 2);
text(0, P.tS+3, 'straight', 'horizontalalignment', 'center', 'fontsize', 7.5);
text(90, P.tL+3, 'left', 'horizontalalignment', 'center', 'fontsize', 7.5);
text(-90, P.tR+3, 'right', 'horizontalalignment', 'center', 'fontsize', 7.5);
text(-165, P.tU-4, 'U', 'horizontalalignment', 'center', 'fontsize', 7.5);
text(165, P.tU-4, 'U', 'horizontalalignment', 'center', 'fontsize', 7.5);
xlim([-180 180]); ylim([-3 46]); set(gca, 'xtick', -180:60:180);
xlabel('Signed deflection angle \theta (deg), + = left'); ylabel('Turn delay (s)');
title('(a) Turn-delay function c(\theta)', 'fontweight', 'normal');
subplot(1,2,2); hold on; box on;
mask = S.K > 0; ang = S.Th(mask); kk = S.K(mask);
ang(ang >= 180 - 1e-9) = 179;   % fold +180 (U-turn) into the last bin
edges = -180:15:180;
for c = 1:4
    hc = histc(ang(kk==c), edges); hc = hc(1:end-1);
    bar(edges(1:end-1)+7.5, hc(:), 1, 'facecolor', turnCol(c,:), 'edgecolor', 'w');
end
xlim([-180 180]); set(gca, 'xtick', -180:60:180);
xlabel('\theta (deg)'); ylabel('Ordered transitions');
title(sprintf('(b) All %d ordered transitions of G^*', sum(mask(:))), 'fontweight', 'normal');
legend(turnName, 'location', 'north', 'fontsize', 7); legend boxoff
savefig(f, fdir, 'fig02_turn_model');

%% fig03 repaired multigraph
f = newfig(6.0, 4.8); hold on; axis equal off;
for e = 1:m
    a = W.E(e,1); b = W.E(e,2);
    plot(xy([a b],1), xy([a b],2), '-', 'color', C.grey, 'linewidth', 2);
end
for e = G.dup'
    a = W.E(e,1); b = W.E(e,2); nv = perp(xy(a,:), xy(b,:)) * 9;
    plot(xy([a b],1)+nv(1), xy([a b],2)+nv(2), '-', 'color', C.red, 'linewidth', 2.4);
end
drawnodes(xy, W.name, R.odd, C);
h1 = plot(nan, nan, '-', 'color', C.grey, 'linewidth', 2);
h2 = plot(nan, nan, '-', 'color', C.red, 'linewidth', 2.4);
legend([h1 h2], {'Street (serviced once)', 'Deadhead copy (T-join)'}, 'location', 'northoutside', ...
    'orientation', 'horizontal', 'fontsize', 8); legend boxoff
xlim([-40 560]); ylim([-40 500]);
savefig(f, fdir, 'fig03_repaired_multigraph');

%% fig04 joint repair analysis
JA = R.joint;
f = newfig(6.2, 4.0); hold on; box on;
for k = 1:size(JA,1)
    plot([JA(k,2) JA(k,2)], [JA(k,3) JA(k,4)], '-', 'color', C.light, 'linewidth', 1);
end
hb = plot(JA(:,2), JA(:,4), 's', 'color', C.grey, 'markersize', 6);
ho = plot(JA(:,2), JA(:,3), 'o', 'markerfacecolor', C.blue, 'markeredgecolor', 'w', 'markersize', 8);
hd = plot(JA(R.kmin,2), JA(R.kmin,3), 'o', 'markeredgecolor', C.orange, 'markersize', 13, 'linewidth', 2);
ht = plot(JA(R.kT,2), JA(R.kT,3), 'd', 'markeredgecolor', C.dark, 'markersize', 13, 'linewidth', 1.5);
xl = [min(JA(:,2))-60, max(JA(:,2))+60];
T0 = JA(R.kT,2)/P.vD + JA(R.kT,3);
for dT = [0 60 120 180]
    xx = linspace(xl(1), xl(2), 50);
    plot(xx, T0 + dT - xx/P.vD, ':', 'color', C.grey);
end
done = false(size(JA,1),1);
for k = 1:size(JA,1)
    if done(k), continue; end
    grp = find(abs(JA(:,2)-JA(k,2)) < 60 & JA(:,3) == JA(k,3) & ~done);
    done(grp) = true;
    lab = strjoin(arrayfun(@(g) sprintf('M%d', g), grp', 'UniformOutput', false), ',');
    if numel(grp) > 1 && max(JA(grp,2)) - min(JA(grp,2)) > 20
        text(mean(JA(grp,2)), JA(k,3)+14, lab, 'fontsize', 6.5, 'color', [0.3 0.3 0.3], 'horizontalalignment', 'center');
    else
        text(JA(k,2)+14, JA(k,3)-9, lab, 'fontsize', 6.5, 'color', [0.3 0.3 0.3]);
    end
end
xlim(xl); ylim([min(JA(:,3))-30, max(JA(:,4))+20]);
xlabel('Deadhead length of the repair (m)'); ylabel('Turn delay of the circuit (s)');
legend([ho hb hd ht], {'Turn-optimal circuit', 'Hierholzer circuit', 'Min-distance repair', 'Min-time repair'}, ...
    'location', 'northoutside', 'orientation', 'horizontal', 'fontsize', 7); legend boxoff
savefig(f, fdir, 'fig04_repair_tradeoff');

%% fig05 baseline vs optimal circuit with turn markers
f = newfig(7.2, 3.9);
subplot(1,2,1); drawcircuit(R.base.seq, R.base.cls, S, G, W, C, turnCol);
title(sprintf('(a) Hierholzer circuit: %d s, %d right, %d U', R.base.cost, sum(R.base.cls==3), sum(R.base.cls==4)), 'fontweight', 'normal');
subplot(1,2,2); drawcircuit(R.opt.seq, R.opt.cls, S, G, W, C, turnCol);
title(sprintf('(b) Turn-optimal circuit: %d s, %d right, %d U', R.opt.cost, sum(R.opt.cls==3), sum(R.opt.cls==4)), 'fontweight', 'normal');
hl = zeros(1,3);
for c = 2:4, hl(c-1) = plot(nan, nan, markerOf(c), 'markerfacecolor', turnCol(c,:), 'markeredgecolor', 'k', 'markersize', 6); end
legend(hl, turnName(2:4), 'location', 'southoutside', 'orientation', 'horizontal', 'fontsize', 7); legend boxoff
savefig(f, fdir, 'fig05_circuits_turns');

%% fig06 distribution of turn delay over uniform random circuits
f = newfig(5.8, 3.4); hold on; box on;
ed = (floor(min(Ms.cost)/12)*12 - 12):12:(ceil(max(Ms.cost)/12)*12 + 12);
hc = histc(Ms.cost, ed); hc = hc(1:end-1);
bar(ed(1:end-1)+6, hc(:)/numel(Ms.cost), 1, 'facecolor', C.light, 'edgecolor', 'w');
yl = [0 1.65*max(hc)/numel(Ms.cost)];
l1 = plot([S.LB S.LB], yl, '--', 'color', C.grey, 'linewidth', 1.6);
l2 = plot([R.opt.cost R.opt.cost], yl, '-', 'color', C.blue, 'linewidth', 1.8);
l3 = plot(mean(Ms.cost)*[1 1], yl, ':', 'color', C.dark, 'linewidth', 1.8);
l4 = plot([R.base.cost R.base.cost], yl, '-', 'color', C.orange, 'linewidth', 1.6);
legend([l1 l2 l3 l4], {sprintf('Lower bound (%d s)', S.LB), sprintf('Exact optimum (%d s)', R.opt.cost), ...
    sprintf('Sample mean (%.0f s)', mean(Ms.cost)), sprintf('Hierholzer (%d s)', R.base.cost)}, ...
    'location', 'northeast', 'fontsize', 7.5); legend boxoff
ylim(yl); xlim([S.LB-30, max(Ms.cost)+20]);
xlabel('Turn delay of the circuit (s)'); ylabel('Fraction of circuits');
title(sprintf('%d uniformly random Eulerian circuits of G^*', numel(Ms.cost)), 'fontweight', 'normal');
savefig(f, fdir, 'fig06_turn_distribution');

%% fig07 turn composition
f = newfig(5.6, 3.3); hold on; box on;
cB = histcounts4(R.base.cls); cO = histcounts4(R.opt.cls);
cR = [mean(Ms.nS) mean(Ms.nL) mean(Ms.nR) mean(Ms.nU)];
sR = [std(Ms.nS) std(Ms.nL) std(Ms.nR) std(Ms.nU)];
Yb = [cB; cR; cO]'; bw = 0.26; cols = [C.orange; C.grey; C.blue];
for s = 1:3
    xs = (1:4) + (s-2)*bw;
    bar(xs, Yb(:,s), bw*0.9, 'facecolor', cols(s,:), 'edgecolor', 'none');
    for k = 1:4
        yt = Yb(k,s) + 0.6; if s == 2, yt = Yb(k,s) + sR(k) + 0.6; end
        text(xs(k), yt, sprintf('%g', round(10*Yb(k,s))/10), 'horizontalalignment', 'center', 'fontsize', 6.5);
    end
end
for k = 1:4, plot([k k], cR(k)+[-1 1]*sR(k), '-', 'color', C.dark, 'linewidth', 0.9); end
set(gca, 'xtick', 1:4, 'xticklabel', turnName); xlim([0.4 4.6]); ylim([0 max(Yb(:))+5]);
ylabel('Number of transitions');
hh = zeros(1,3); for s = 1:3, hh(s) = plot(nan, nan, 's', 'markerfacecolor', cols(s,:), 'markeredgecolor', 'none', 'markersize', 8); end
legend(hh, {'Hierholzer', 'Random mean (\pm1 sd)', 'Turn-optimal'}, 'location', 'northeast', 'fontsize', 7.5); legend boxoff
savefig(f, fdir, 'fig07_turn_composition');

%% fig08 trips
f = newfig(7.2, 3.9);
tripCol = [C.blue; C.orange; C.aqua; C.violet; C.yellow; C.red];
subplot(1,2,1); drawtrips(R.base.seq, Sb, S, G, W, C, tripCol);
title(sprintf('(a) Hierholzer: %d trips, %.1f min', Sb.ntrips, Sb.total/60), 'fontweight', 'normal');
subplot(1,2,2); drawtrips(R.opt.seq, So, S, G, W, C, tripCol);
title(sprintf('(b) Turn-optimal: %d trips, %.1f min', So.ntrips, So.total/60), 'fontweight', 'normal');
hh = zeros(1, So.ntrips); lab = cell(1, So.ntrips);
for t = 1:So.ntrips
    hh(t) = plot(nan, nan, '-', 'color', tripCol(t,:), 'linewidth', 2.4);
    lab{t} = sprintf('Trip %d', t);
end
legend(hh, lab, 'location', 'southoutside', 'orientation', 'horizontal', 'fontsize', 7); legend boxoff
savefig(f, fdir, 'fig08_trips');

%% fig09 shift-time breakdown
rs = R.splitRand;
f = newfig(7.0, 3.0);
subplot(1,2,1); hold on; box on;
Mx = [Sb.tService Sb.tDeadIn Sb.tAccess Sb.tTurn Sb.tUnload;
      mean(rs(:,2:6));
      So.tService So.tDeadIn So.tAccess So.tTurn So.tUnload] / 60;
segCol = [C.blue; C.orange; C.yellow; C.red; C.grey];
for r = 1:3
    x0 = 0;
    for k = 1:5
        patch([x0 x0+Mx(r,k) x0+Mx(r,k) x0], [r-0.32 r-0.32 r+0.32 r+0.32], segCol(k,:), 'edgecolor', 'w');
        x0 = x0 + Mx(r,k);
    end
    text(x0+2, r, sprintf('%.1f', x0), 'fontsize', 7.5);
end
set(gca, 'ytick', 1:3, 'yticklabel', {'Hierholzer', 'Random mean', 'Turn-optimal'}, 'ydir', 'reverse');
xlim([0 1.18*max(sum(Mx,2))]); ylim([0.4 3.6]); xlabel('Shift time (min)');
hh = zeros(1,5); for k = 1:5, hh(k) = patch(nan, nan, segCol(k,:)); end
legend(hh, {'Service', 'In-ward deadhead', 'Access', 'Turns', 'Unloading'}, 'location', 'southoutside', ...
    'orientation', 'horizontal', 'fontsize', 6.5); legend boxoff
title('(a) Composition of the daily shift', 'fontweight', 'normal');
subplot(1,2,2); hold on; box on;
ed = floor(min(rs(:,1))/60):0.25:ceil(max(rs(:,1))/60)+0.25;
hc = histc(rs(:,1)/60, ed); hc = hc(1:end-1);
bar(ed(1:end-1)+0.125, hc(:), 1, 'facecolor', C.light, 'edgecolor', 'w');
yl = [0 1.2*max(hc)];
vl(So.total/60, yl, C.blue, '-', 'optimum');
vl(Sb.total/60, yl, C.orange, '-', 'Hierholzer');
ylim(yl); xlabel('Shift time (min)'); ylabel('Random circuits');
title(sprintf('(b) Shift time of %d random circuits', size(rs,1)), 'fontweight', 'normal');
savefig(f, fdir, 'fig09_shift_time');

%% fig10 sensitivity
f = newfig(7.0, 2.9);
sv = R.sens;
subplot(1,2,1); heat(sv.red, sv.tR, sv.tL, '%.0f', [1 1 1; C.blue]);
xlabel('Right-turn delay t_R (s)'); ylabel('Left-turn delay t_L (s)');
title('(a) Turn-delay reduction vs Hierholzer (%)', 'fontweight', 'normal');
subplot(1,2,2); heat(sv.nR, sv.tR, sv.tL, '%d', [1 1 1; C.orange]);
xlabel('Right-turn delay t_R (s)'); ylabel('Left-turn delay t_L (s)');
title(sprintf('(b) Right turns in the optimum (Hierholzer: %d)', sv.nRbase(1)), 'fontweight', 'normal');
savefig(f, fdir, 'fig10_penalty_sensitivity');

%% fig11 capacity sensitivity
cs = R.capsens;
f = newfig(7.0, 2.9);
subplot(1,2,1); hold on; box on;
plot(cs(:,1), cs(:,2)/60, '-s', 'color', C.orange, 'markersize', 4, 'markerfacecolor', C.orange, 'linewidth', 1.5);
plot(cs(:,1), cs(:,3)/60, '-o', 'color', C.blue, 'markersize', 4, 'markerfacecolor', C.blue, 'linewidth', 1.5);
ix = find(cs(:,8) > 0 | cs(:,6) ~= cs(:,3));
plot(cs(ix,1), cs(ix,6)/60, 'd', 'color', C.dark, 'markersize', 7, 'linewidth', 1.5);
xlabel('Truck payload Q (kg)'); ylabel('Shift time (min)');
legend({'Hierholzer', 'Turn-optimal', 'Split-aware LS'}, 'location', 'northeast', 'fontsize', 7.5); legend boxoff
title('(a) Shift time', 'fontweight', 'normal');
subplot(1,2,2); hold on; box on;
bar(cs(:,1), (cs(:,2)-cs(:,3))/60, 0.8, 'facecolor', C.light, 'edgecolor', 'none');
bar(cs(:,1), (cs(:,2)-cs(:,6))/60, 0.4, 'facecolor', C.blue, 'edgecolor', 'none');
plot(cs([1 end],1), [0 0], '-', 'color', C.dark);
ch = find(diff(cs(:,7)) ~= 0);
yl = ylim;
for k = [0; ch]'
    text(cs(k+1,1), yl(2)*0.92, sprintf('%d', cs(k+1,7)), 'fontsize', 7, 'color', C.dark);
end
xlabel('Truck payload Q (kg)'); ylabel('Saving (min)');
legend({'Turn-optimal', 'With split-aware LS'}, 'location', 'south', 'fontsize', 7); legend boxoff
title('(b) Saving over Hierholzer (labels: trips)', 'fontweight', 'normal');
savefig(f, fdir, 'fig11_capacity_sensitivity');

%% fig12 benchmark
if doBench && isfield(R, 'bench')
    BM = R.bench; sz = R.benchSizes; ns = size(sz,1);
    f = newfig(7.0, 3.0);
    subplot(1,2,1); hold on; box on;
    gh = 100*(BM(:,10)./BM(:,11) - 1); gb = 100*(BM(:,9)./BM(:,11) - 1);
    jit = 0.12*((BM(:,2)-5.5)/5);
    plot(BM(:,1)-0.15+jit, gb, 's', 'color', C.orange, 'markersize', 4);
    plot(BM(:,1)+0.15+jit, gh, 'o', 'color', C.blue, 'markersize', 4);
    mb = zeros(ns,1); mh = zeros(ns,1);
    for s = 1:ns, mb(s) = mean(gb(BM(:,1)==s)); mh(s) = mean(gh(BM(:,1)==s)); end
    h1 = plot((1:ns)-0.15, mb, '-', 'color', C.orange, 'linewidth', 1.8);
    h2 = plot((1:ns)+0.15, mh, '-', 'color', C.blue, 'linewidth', 1.8);
    lab = arrayfun(@(s) sprintf('%dx%d', sz(s,1), sz(s,2)), 1:ns, 'UniformOutput', false);
    set(gca, 'xtick', 1:ns, 'xticklabel', lab); xlim([0.5 ns+0.5]);
    xlabel('Grid size'); ylabel('Excess turn delay over exact (%)');
    legend([h1 h2], {'Hierholzer', 'Merge + local search'}, 'location', 'northwest', 'fontsize', 7.5); legend boxoff
    title('(a) Optimality gap (10 instances per size)', 'fontweight', 'normal');
    subplot(1,2,2); hold on; box on;
    nn = zeros(ns,1); ti = zeros(ns,1); th_ = zeros(ns,1); tmx = zeros(ns,1);
    for s = 1:ns
        b = BM(BM(:,1)==s,:);
        nn(s) = mean(b(:,5)); ti(s) = median(b(:,12)); th_(s) = median(b(:,13)); tmx(s) = max(b(:,12));
    end
    loglog(nn, ti, '-o', 'color', C.dark, 'markerfacecolor', C.dark, 'markersize', 5, 'linewidth', 1.5);
    loglog(nn, tmx, ':', 'color', C.dark, 'linewidth', 1);
    loglog(nn, th_, '-s', 'color', C.blue, 'markerfacecolor', C.blue, 'markersize', 5, 'linewidth', 1.5);
    set(gca, 'xscale', 'log', 'yscale', 'log');
    xlabel('Edges of the repaired multigraph m^*'); ylabel('Time (s)');
    legend({'Exact 0-1 programme (median)', 'Exact (max)', 'Heuristic (median)'}, 'location', 'northwest', 'fontsize', 7.5); legend boxoff
    title('(b) Running time', 'fontweight', 'normal');
    savefig(f, fdir, 'fig12_benchmark');
end
end

%% =================================================================== helpers
function f = newfig(w, h)
f = figure('visible', 'off', 'color', 'w', 'units', 'inches', 'position', [0 0 w h]);
set(f, 'paperunits', 'inches', 'papersize', [w h], 'paperposition', [0 0 w h]);
end

function savefig(f, fdir, name)
print(f, fullfile(fdir, [name '.pdf']), '-dpdf');
print(f, fullfile(fdir, [name '.png']), '-dpng', '-r200');
close(f);
end

function nv = perp(p, q)
d = q - p; d = d / norm(d); nv = [-d(2) d(1)];
end

function drawnodes(xy, names, odd, C)
n = size(xy,1);
for v = 1:n
    if any(odd == v), fc = C.red; else, fc = 'w'; end
    plot(xy(v,1), xy(v,2), 'o', 'markersize', 7, 'markerfacecolor', fc, 'markeredgecolor', 'k');
    text(xy(v,1)-14, xy(v,2)+17, names{v}, 'fontsize', 7.5, 'fontweight', 'bold');
end
end

function off = edgeoffset(G, e)
% perpendicular offset (m) so that the two copies of a duplicated street are both visible
if ~G.svc(e), off = 9;
elseif any(G.dup == G.orig(e)), off = -9;
else, off = 0; end
end

function drawcircuit(seq, cls, S, G, W, C, turnCol)
hold on; axis equal off;
xy = W.xy; L = numel(seq);
ramp = [linspace(0.80, 0.05, L)' linspace(0.88, 0.25, L)' linspace(0.98, 0.55, L)'];
for k = 1:L
    h = seq(k); e = ceil(h/2);
    a = S.hv(h); b = S.hv(S.oth(h));
    nv = perp(xy(min(a,b),:), xy(max(a,b),:)) * edgeoffset(G, e);
    p = xy(a,:) + nv; q = xy(b,:) + nv;
    plot([p(1) q(1)], [p(2) q(2)], '-', 'color', ramp(k,:), 'linewidth', 2);
    mid = 0.5*(p+q); d = (q-p)/norm(q-p);
    plot([mid(1)-9*d(1)-6*d(2), mid(1)+9*d(1), mid(1)-9*d(1)+6*d(2)], ...
         [mid(2)-9*d(2)+6*d(1), mid(2)+9*d(2), mid(2)-9*d(2)-6*d(1)], '-', 'color', ramp(k,:), 'linewidth', 1.2);
end
for v = 1:size(xy,1)
    plot(xy(v,1), xy(v,2), 'o', 'markersize', 4, 'markerfacecolor', 'w', 'markeredgecolor', C.dark);
end
% turn markers, offset into the turn
for k = 1:L
    c = cls(k);
    if c == 1, continue; end
    arr = S.oth(seq(k)); nxt = seq(mod(k,L)+1);
    v = S.hv(arr);
    din = xy(v,:) - xy(S.hv(S.oth(arr)),:); din = din/norm(din);
    dout = xy(S.hv(S.oth(nxt)),:) - xy(v,:); dout = dout/norm(dout);
    o = dout - din; if norm(o) < 1e-6, o = [0 1]; end
    o = o/norm(o) * (16 + 8*mod(k,2));
    plot(xy(v,1)+o(1), xy(v,2)+o(2), markerOf(c), 'markerfacecolor', turnCol(c,:), ...
        'markeredgecolor', 'k', 'markersize', 5);
end
s0 = S.hv(seq(1));
text(xy(s0,1)-40, xy(s0,2)-28, 'start', 'fontsize', 7, 'color', C.dark);
xlim([-50 560]); ylim([-50 500]);
end

function mk = markerOf(c)
mk = {'.', 'o', '^', 's'}; mk = mk{c};
end

function drawtrips(seq, Sx, S, G, W, C, tripCol)
hold on; axis equal off;
xy = W.xy;
for e = 1:size(W.E,1)
    plot(xy(W.E(e,:),1), xy(W.E(e,:),2), '-', 'color', [0.92 0.92 0.91], 'linewidth', 4);
end
for t = 1:Sx.ntrips
    pos = Sx.trip(t).pos;
    for k = pos
        h = seq(k); e = ceil(h/2);
        a = S.hv(h); b = S.hv(S.oth(h));
        nv = perp(xy(min(a,b),:), xy(max(a,b),:)) * edgeoffset(G, e);
        if G.svc(e), ls = '-'; lw = 2.4; else, ls = '--'; lw = 1.4; end
        plot(xy([a b],1)+nv(1), xy([a b],2)+nv(2), ls, 'color', tripCol(t,:), 'linewidth', lw);
    end
    a = S.hv(seq(pos(1))); b = S.hv(S.oth(seq(pos(end))));
    plot(xy(a,1), xy(a,2), 'o', 'markersize', 8, 'markeredgecolor', tripCol(t,:), 'linewidth', 2);
    plot(xy(b,1), xy(b,2), 'x', 'markersize', 9, 'color', tripCol(t,:), 'linewidth', 2);
end
for v = 1:size(xy,1)
    plot(xy(v,1), xy(v,2), '.', 'color', C.dark, 'markersize', 8);
end
xlim([-50 560]); ylim([-50 500]);
end

function vl(x, yl, col, ls, lab)
plot([x x], yl, ls, 'color', col, 'linewidth', 1.6);
text(x, yl(2)*0.97, [' ' lab], 'color', [0.15 0.15 0.15], 'fontsize', 7, 'rotation', 90, ...
    'horizontalalignment', 'right', 'verticalalignment', 'bottom');
end

function c = histcounts4(cls)
c = [sum(cls==1) sum(cls==2) sum(cls==3) sum(cls==4)];
end

function heat(M, xv, yv, fmt, cl)
mx = max(M(:)); mn = min(M(:)); if mx == mn, mx = mn + 1; end
hold on;
for i = 1:size(M,1)
    for j = 1:size(M,2)
        t = (M(i,j) - mn)/(mx - mn);
        col = cl(1,:)*(1-t) + cl(2,:)*t;
        patch([j-0.5 j+0.5 j+0.5 j-0.5], [i-0.5 i-0.5 i+0.5 i+0.5], col, 'edgecolor', 'w');
        if t > 0.6, tc = 'w'; else, tc = 'k'; end
        text(j, i, sprintf(fmt, M(i,j)), 'horizontalalignment', 'center', 'fontsize', 7.5, 'color', tc);
    end
end
set(gca, 'xtick', 1:numel(xv), 'xticklabel', arrayfun(@num2str, xv, 'UniformOutput', false), ...
    'ytick', 1:numel(yv), 'yticklabel', arrayfun(@num2str, yv, 'UniformOutput', false));
xlim([0.5 numel(xv)+0.5]); ylim([0.5 numel(yv)+0.5]); box on;
end
