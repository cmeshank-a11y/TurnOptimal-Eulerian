function G = eulerize_closed(n, E, len, pairs, Nx)
%EULERIZE_CLOSED  Duplicate the streets on the shortest paths that join the
%   matched odd vertices (a minimum T-join when pairs is a minimum-weight
%   matching). Returns the Eulerian multigraph G with fields
%   E (edges), len, svc (true = service copy), orig (street index).
m = size(E,1);
key = E(:,1)*(n+1) + E(:,2);
dup = [];
for r = 1:size(pairs,1)
    p = sp_path(Nx, pairs(r,1), pairs(r,2));
    for s = 1:numel(p)-1
        u = min(p(s),p(s+1)); v = max(p(s),p(s+1));
        cand = find(key == u*(n+1)+v);
        [~, b] = min(len(cand));
        dup(end+1) = cand(b); %#ok<AGROW>
    end
end
% a street duplicated twice is restored to parity by removing both copies
dup = sort(dup);
keep = true(size(dup));
for t = 1:numel(dup)-1
    if keep(t) && keep(t+1) && dup(t) == dup(t+1)
        keep(t) = false; keep(t+1) = false;
    end
end
dup = dup(keep);
G.n = n;
G.E = [E; E(dup,:)];
G.len = [len(:); len(dup(:))];
G.svc = [true(m,1); false(numel(dup),1)];
G.orig = [(1:m)'; dup(:)];
G.dup = dup(:);
G.m = size(G.E,1);
end
