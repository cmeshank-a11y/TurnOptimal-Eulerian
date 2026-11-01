function [cost, seq, ok, Fw, Rv] = ts_eval(S, partner)
%TS_EVAL  Exact orientation-optimal turn cost of a transition system.
%   ok = true when it is a single Eulerian circuit. seq is the sequence of
%   departing half-edges in the cheaper orientation (cost = min(F,R)).
[L, Fw, Rv, seq] = ts_trace(S, partner, 1);
ok = (L == S.m);
if ~ok, cost = inf; return; end
if Rv < Fw
    cost = Rv;
    seq = S.oth(flipud(seq));
else
    cost = Fw;
end
end
