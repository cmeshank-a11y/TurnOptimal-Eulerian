function [partner, eseq, vseq] = hierholzer_partner(S, G, start)
%HIERHOLZER_PARTNER  Classical Hierholzer circuit (incidence lists in
%   increasing edge order, i.e. what a distance-only planner returns) and
%   its transition system. eseq(k) is the k-th edge, vseq(k) its tail.
n = G.n; m = G.m;
adj = cell(n,1);
for e = 1:m
    adj{G.E(e,1)}(end+1) = e; adj{G.E(e,2)}(end+1) = e;
end
ptr = ones(n,1); used = false(m,1);
vs = start; es = []; cv = []; ce = [];
while ~isempty(vs)
    v = vs(end);
    while ptr(v) <= numel(adj{v}) && used(adj{v}(ptr(v)))
        ptr(v) = ptr(v) + 1;
    end
    if ptr(v) <= numel(adj{v})
        e = adj{v}(ptr(v)); used(e) = true;
        u = G.E(e,1) + G.E(e,2) - v;
        vs(end+1) = u; es(end+1) = e; %#ok<AGROW>
    else
        cv(end+1) = v; vs(end) = []; %#ok<AGROW>
        if ~isempty(es), ce(end+1) = es(end); es(end) = []; end %#ok<AGROW>
    end
end
vseq = fliplr(cv); eseq = fliplr(ce);        % vseq(k) -> vseq(k+1) via eseq(k)
vseq = vseq(1:m);
partner = zeros(S.H,1);
for k = 1:m
    kn = mod(k, m) + 1;
    v = vseq(kn);                            % vertex between eseq(k) and eseq(kn)
    a = he(G, eseq(k), v); b = he(G, eseq(kn), v);
    partner(a) = b; partner(b) = a;
end
end
function h = he(G, e, v)
if G.E(e,1) == v, h = 2*e - 1; else, h = 2*e; end
end
