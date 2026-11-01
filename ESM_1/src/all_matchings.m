function M = all_matchings(T)
%ALL_MATCHINGS  Every perfect matching of the list T, as a cell array of
%   (numel(T)/2)-by-2 arrays. There are (numel(T)-1)!! of them.
if isempty(T), M = {zeros(0,2)}; return; end
M = {};
a = T(1);
for j = 2:numel(T)
    rest = T([2:j-1, j+1:end]);
    S = all_matchings(rest);
    for s = 1:numel(S)
        M{end+1} = [a T(j); S{s}]; %#ok<AGROW>
    end
end
end
