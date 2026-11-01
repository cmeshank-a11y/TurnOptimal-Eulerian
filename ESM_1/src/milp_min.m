function [x, fval, ok] = milp_min(c, Aeq, beq, Age, bge, tlim)
%MILP_MIN  Solve  min c'x  s.t.  Aeq x = beq,  Age x >= bge,  x binary.
%   Uses intlinprog (MATLAB, Optimization Toolbox) when available and GLPK
%   (built into GNU Octave) otherwise. tlim: time limit in seconds (default inf);
%   ok = false if the limit is hit before optimality is proven.
if nargin < 6, tlim = inf; end
n = numel(c); lb = zeros(n,1); ub = ones(n,1);
isOct = exist('OCTAVE_VERSION', 'builtin') > 0;
if ~isOct && exist('intlinprog', 'file') == 2
    opts = optimoptions('intlinprog', 'Display', 'off');
    if isfinite(tlim), opts = optimoptions(opts, 'MaxTime', tlim); end
    [x, fval, flag] = intlinprog(c(:), 1:n, -Age, -bge, Aeq, beq, lb, ub, opts);
    ok = flag == 1;
else
    A = [Aeq; Age]; b = [beq; bge];
    ctype = [repmat('S', 1, size(Aeq,1)), repmat('L', 1, size(Age,1))];
    param.msglev = 0;
    if isfinite(tlim), param.tmlim = round(1000*tlim); end
    [x, fval, err, extra] = glpk(c(:), A, b, lb, ub, ctype, repmat('I',1,n), 1, param);
    ok = (err == 0) && (extra.status == 5);
end
x = round(x);
end
