function [partner, Sx, nmoves] = split_ls(S, G, partner, W, D, P, Q, maxPass)
%SPLIT_LS  Split-aware local search: starting from a circuit, re-pair the
%   transitions of one vertex at a time and keep the change when the result
%   is still one Eulerian circuit and its optimally split shift time
%   (SPLIT_TOUR over all rotations and orientations) decreases.
if nargin < 8, maxPass = 5; end
Sx = best_split(S, G, partner, W, D, P, Q);
vary = find(S.np > 1)';
nmoves = 0;
for pass = 1:maxPass
    improved = false;
    for v = vary
        for p = 1:S.np(v)
            trial = ts_apply(partner, S.pairs{v}{p});
            [L, ~, ~] = ts_trace(S, trial, 1);
            if L ~= S.m, continue; end
            St = best_split(S, G, trial, W, D, P, Q);
            if St.total < Sx.total - 1e-6
                partner = trial; Sx = St; improved = true; nmoves = nmoves + 1;
            end
        end
    end
    if ~improved, break; end
end
end

function Sx = best_split(S, G, partner, W, D, P, Q)
% both orientations of the circuit are split; the faster one is kept
[~, ~, ~, seq] = ts_trace(S, partner, 1);
A = split_tour(S, G, seq, W, D, P, Q);
B = split_tour(S, G, S.oth(flipud(seq)), W, D, P, Q);
if B.total < A.total, Sx = B; Sx.seq = S.oth(flipud(seq)); else, Sx = A; Sx.seq = seq; end
end
