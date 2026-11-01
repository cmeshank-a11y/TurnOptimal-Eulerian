function S = ts_prepare(G, xy, P)
%TS_PREPARE  Half-edge structure, turn-cost matrix and all transition
%   pairings of the Eulerian multigraph G with vertex coordinates xy.
%   Half-edge 2e-1 is the end of edge e at G.E(e,1), half-edge 2e the end
%   at G.E(e,2). A transition (a,b) at v means: arrive at v along the edge
%   of half-edge a, depart along the edge of half-edge b.
%   P: turn penalties in seconds (tS,tL,tR,tU) and thresholds (degS,degU).
n = G.n; m = G.m; H = 2*m;
hv = zeros(H,1); hv(1:2:H) = G.E(:,1); hv(2:2:H) = G.E(:,2);
oth = zeros(H,1); oth(1:2:H) = 2:2:H; oth(2:2:H) = 1:2:H;
inc = cell(n,1);
for h = 1:H, inc{hv(h)}(end+1) = h; end
C = zeros(H); K = zeros(H); Th = zeros(H);
pen = [P.tS P.tL P.tR P.tU];
for v = 1:n
    for a = inc{v}
        din = xy(v,:) - xy(hv(oth(a)),:);
        for b = inc{v}
            if a == b, continue; end
            dout = xy(hv(oth(b)),:) - xy(v,:);
            th = atan2(din(1)*dout(2) - din(2)*dout(1), din*dout') * 180/pi;
            if abs(th) <= P.degS
                c = 1;                 % straight
            elseif abs(th) > P.degU
                c = 4;                 % U-turn
            elseif th > 0
                c = 2;                 % left turn (counter-clockwise)
            else
                c = 3;                 % right turn: crosses opposing flow (LHT)
            end
            Th(a,b) = th; K(a,b) = c; C(a,b) = pen(c);
        end
    end
end
Cmin = min(C, C');                     % cheaper of the two directions
pairs = cell(n,1); mc = cell(n,1); np = zeros(n,1);
for v = 1:n
    L = all_matchings(inc{v});
    c = zeros(numel(L),1);
    for p = 1:numel(L)
        Q = L{p};
        c(p) = sum(Cmin(sub2ind([H H], Q(:,1), Q(:,2))));
    end
    [c, o] = sort(c);                  % ascending: B&B relies on this order
    pairs{v} = L(o); mc{v} = c; np(v) = numel(L);
end
S.n = n; S.m = m; S.H = H; S.hv = hv; S.oth = oth; S.inc = inc;
S.C = C; S.K = K; S.Th = Th; S.Cmin = Cmin;
S.pairs = pairs; S.mc = mc; S.np = np;
S.LB = sum(cellfun(@(x) x(1), mc));    % Proposition: local lower bound
S.logNTS = sum(log10(np));             % log10 of the number of transition systems
end
