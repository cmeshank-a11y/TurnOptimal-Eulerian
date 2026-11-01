function [L, F, R, seq] = ts_trace(S, partner, h0)
%TS_TRACE  Follow the closed trail of the transition system 'partner' that
%   departs along half-edge h0. L = number of edges on it (L == S.m iff the
%   transition system is a single Eulerian circuit), F and R = turn cost in
%   the forward and reverse orientation, seq = departing half-edges.
if nargin < 3, h0 = 1; end
dep = h0; L = 0; F = 0; R = 0;
want = nargout > 3;
if want, seq = zeros(S.m,1); end
while true
    L = L + 1;
    if want, seq(L) = dep; end
    arr = S.oth(dep);
    nxt = partner(arr);
    F = F + S.C(arr, nxt);
    R = R + S.C(nxt, arr);
    dep = nxt;
    if dep == h0 || L > S.m, break; end
end
if want, seq = seq(1:min(L,S.m)); end
end
