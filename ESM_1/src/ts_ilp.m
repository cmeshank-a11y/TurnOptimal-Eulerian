function X = ts_ilp(S, G, maxRounds, tlim)
%TS_ILP  Exact minimum turn-cost Eulerian circuit as a 0-1 programme over
%   edge orientations y_e and directed transitions x_(a->b), with lazy
%   connectivity cuts (Proposition 4.5):
%     sum_b x_(h->b) = [edge of h arrives at v],
%     sum_a x_(a->h) = [edge of h departs from v]   for every half-edge h,
%     sum over transitions leaving edge set U  >= 1  for every closed trail U
%   found in an intermediate solution that is not the whole edge set.
if nargin < 3 || isempty(maxRounds), maxRounds = 500; end
if nargin < 4, tlim = inf; end   % total time budget (s)
t0 = tic;
m = S.m; H = S.H;
% enumerate directed transitions
ta = []; tb = [];
for v = 1:S.n
    hs = S.inc{v};
    for a = hs, for b = hs
        if a ~= b, ta(end+1,1) = a; tb(end+1,1) = b; end %#ok<AGROW>
    end, end
end
T = numel(ta);
nv = m + T;                                 % [y ; x]
c = [zeros(m,1); S.C(sub2ind([H H], ta, tb))];
Aeq = zeros(2*H, nv); beq = zeros(2*H, 1);
for h = 1:H
    e = ceil(h/2); odd = mod(h,2) == 1;
    r1 = 2*h - 1; r2 = 2*h;
    Aeq(r1, m + find(ta == h)) = 1;         % leave h's vertex after arriving via h
    Aeq(r2, m + find(tb == h)) = 1;         % arrive somewhere, depart via h
    if odd                                  % h at tail E(e,1): arrive iff y_e = 0
        Aeq(r1, e) = 1;  beq(r1) = 1;
        Aeq(r2, e) = -1; beq(r2) = 0;
    else                                    % h at head E(e,2): arrive iff y_e = 1
        Aeq(r1, e) = -1; beq(r1) = 0;
        Aeq(r2, e) = 1;  beq(r2) = 1;
    end
end
Age = zeros(0, nv); bge = zeros(0,1);
ea = ceil(ta/2); eb = ceil(tb/2);
X.rounds = 0; X.cuts = 0; X.ok = false; X.cost = NaN; X.partner = [];
for round_ = 1:maxRounds
    left = tlim - toc(t0);
    if left <= 0, break; end
    [x, fval, ok] = milp_min(c, Aeq, beq, Age, bge, left);
    if ~ok, break; end
    X.rounds = round_;
    sel = find(x(m+1:end) > 0.5);
    partner = zeros(H,1);
    partner(ta(sel)) = tb(sel); partner(tb(sel)) = ta(sel);
    [comp, nc] = ts_components(S, partner);
    if nc == 1
        X.ok = true; X.partner = partner; X.cost = fval;
        X.y = x(1:m); X.x = x(m+1:end);
        break;
    end
    for k = 1:nc                            % one cut per closed trail
        inU = comp == k;
        row = zeros(1, nv);
        row(m + find(inU(ea) & ~inU(eb))) = 1;
        Age(end+1,:) = row; bge(end+1,1) = 1; %#ok<AGROW>
        X.cuts = X.cuts + 1;
    end
end
X.nvar = nv; X.ntrans = T; X.time = toc(t0);
end
