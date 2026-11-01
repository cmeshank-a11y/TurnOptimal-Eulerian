function u = wh_rand(arg, k)
%WH_RAND  Wichmann-Hill (AS 183) uniform generator, exact in double precision.
%   wh_rand('seed', k)  re-seeds the generator from the integer k >= 0.
%   u = wh_rand(n)      returns an n-by-1 vector of U(0,1) numbers.
%   The generator is coded here (not MATLAB's rand) so that every random
%   quantity in the study is bit-identical in MATLAB and GNU Octave.
persistent s
if isempty(s), s = [1 2 3]; end
if ischar(arg)
    s = [1 + mod(k, 30000), 1 + mod(3*k + 7, 30000), 1 + mod(7*k + 11, 30000)];
    u = [];
    return
end
n = arg;
u = zeros(n, 1);
for i = 1:n
    s(1) = mod(171*s(1), 30269);
    s(2) = mod(172*s(2), 30307);
    s(3) = mod(170*s(3), 30323);
    u(i) = mod(s(1)/30269 + s(2)/30307 + s(3)/30323, 1);
end
end
