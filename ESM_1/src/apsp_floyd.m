function [D, Nx] = apsp_floyd(n, E, w)
%APSP_FLOYD  All-pairs shortest paths (Floyd-Warshall) with successor matrix.
%   Nx(i,j) is the vertex after i on a shortest i-j path.
D = inf(n);
D(1:n+1:end) = 0;
Nx = zeros(n);
for i = 1:n, Nx(i,i) = i; end
for k = 1:size(E,1)
    u = E(k,1); v = E(k,2);
    if w(k) < D(u,v)
        D(u,v) = w(k); D(v,u) = w(k);
        Nx(u,v) = v;   Nx(v,u) = u;
    end
end
for k = 1:n
    M = D(:,k) + D(k,:);
    mask = M < D - 1e-9;
    Nk = repmat(Nx(:,k), 1, n);
    D(mask) = M(mask);
    Nx(mask) = Nk(mask);
end
end
