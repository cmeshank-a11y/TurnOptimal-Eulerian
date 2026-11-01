function c = ts_mincost(S, partner)
%TS_MINCOST  Sum over transitions of the cheaper-direction turn cost
%   (the relaxation used in the lower bound and the merge heuristic).
h = (1:S.H)';
sel = h < partner(h);
c = sum(S.Cmin(sub2ind([S.H S.H], h(sel), partner(h(sel)))));
end
