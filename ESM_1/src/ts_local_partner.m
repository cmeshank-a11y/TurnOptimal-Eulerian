function partner = ts_local_partner(S, idx)
%TS_LOCAL_PARTNER  Transition system choosing pairing idx(v) at each vertex
%   (idx = ones gives the vertex-wise cheapest pairing).
partner = zeros(S.H,1);
for v = 1:S.n
    partner = ts_apply(partner, S.pairs{v}{idx(v)});
end
end
