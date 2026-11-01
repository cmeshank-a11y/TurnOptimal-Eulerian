function [comp, nc] = ts_components(S, partner)
%TS_COMPONENTS  Label each edge by the closed trail containing it.
comp = zeros(S.m,1); nc = 0;
for e = 1:S.m
    if comp(e), continue; end
    nc = nc + 1;
    dep = 2*e - 1;
    while true
        comp(ceil(dep/2)) = nc;
        dep = partner(S.oth(dep));
        if dep == 2*e - 1, break; end
    end
end
end
