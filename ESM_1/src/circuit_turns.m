function [cls, cst, ang] = circuit_turns(S, seq)
%CIRCUIT_TURNS  Turn class (1 S, 2 L, 3 R, 4 U), delay and signed angle of
%   every transition of the oriented circuit given by departing half-edges.
m = numel(seq);
cls = zeros(m,1); cst = zeros(m,1); ang = zeros(m,1);
for k = 1:m
    arr = S.oth(seq(k)); nxt = seq(mod(k,m)+1);
    cls(k) = S.K(arr,nxt); cst(k) = S.C(arr,nxt); ang(k) = S.Th(arr,nxt);
end
end
