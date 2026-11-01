function [pairs, cost] = min_matching_dp(T, D)
%MIN_MATCHING_DP  Exact minimum-weight perfect matching of the vertex set T
%   under the metric D, by dynamic programming over subsets (O(2^k k)).
k = numel(T);
if k == 0, pairs = zeros(0,2); cost = 0; return; end
N = 2^k;
f = inf(N,1); f(1) = 0;
prv = zeros(N,1); pi_ = zeros(N,1); pj_ = zeros(N,1);
bits = 2.^(0:k-1);
for mask = 0:N-1
    if ~isfinite(f(mask+1)), continue; end
    i = find(bitand(mask, bits) == 0, 1);
    if isempty(i), continue; end
    for j = i+1:k
        if bitand(mask, bits(j)), continue; end
        nm = mask + bits(i) + bits(j);
        c = f(mask+1) + D(T(i), T(j));
        if c < f(nm+1) - 1e-9
            f(nm+1) = c; prv(nm+1) = mask; pi_(nm+1) = i; pj_(nm+1) = j;
        end
    end
end
cost = f(N);
pairs = zeros(k/2, 2);
mask = N - 1; r = 0;
while mask > 0
    r = r + 1;
    pairs(r,:) = [T(pi_(mask+1)), T(pj_(mask+1))];
    mask = prv(mask+1);
end
end
