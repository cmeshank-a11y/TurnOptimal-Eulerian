function p = sp_path(Nx, i, j)
%SP_PATH  Vertex sequence of the shortest i-j path encoded by Nx.
p = i;
while p(end) ~= j
    p(end+1) = Nx(p(end), j); %#ok<AGROW>
end
end
