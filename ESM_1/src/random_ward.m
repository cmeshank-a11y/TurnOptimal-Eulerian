function W = random_ward(r, c, pdiag, seed)
%RANDOM_WARD  Random planar street grid of r rows and c columns: jittered
%   150 m blocks, every block crossed by one diagonal lane with probability
%   pdiag, random street classes. Used for the benchmark of Section 7.6.
wh_rand('seed', seed);
n = r*c;
[I, J] = meshgrid(1:c, 1:r);              % J = row, I = column
I = I'; J = J';
xy = [150*(I(:)-1), 150*(J(:)-1)] + 40*(reshape(wh_rand(2*n), n, 2) - 0.5);
id = @(i,j) (j-1)*c + i;
E = zeros(0,2);
for j = 1:r, for i = 1:c-1, E(end+1,:) = [id(i,j) id(i+1,j)]; end, end %#ok<AGROW>
for j = 1:r-1, for i = 1:c, E(end+1,:) = [id(i,j) id(i,j+1)]; end, end %#ok<AGROW>
for j = 1:r-1
    for i = 1:c-1
        u = wh_rand(2);
        if u(1) < pdiag
            if u(2) < 0.5, E(end+1,:) = [id(i,j) id(i+1,j+1)]; %#ok<AGROW>
            else,          E(end+1,:) = [id(i+1,j) id(i,j+1)]; end %#ok<AGROW>
        end
    end
end
E = sort(E, 2);
W.xy = xy; W.E = E;
W.len = round(sqrt(sum((xy(E(:,1),:) - xy(E(:,2),:)).^2, 2)));
W.class = 1 + floor(3*wh_rand(size(E,1)));
W.rho = [0.55 0.85 1.30];
W.demand = W.len .* W.rho(W.class)';
W.gate = [1 c]; W.gateLen = [1200 1200];
W.name = arrayfun(@(k) sprintf('v%d',k), 1:n, 'UniformOutput', false);
end
