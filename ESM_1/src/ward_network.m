function W = ward_network()
%WARD_NETWORK  Illustrative 16-intersection, 27-street kerbside-collection ward.
%   Coordinates are in metres (x east, y north). Street lengths are the
%   Euclidean distances between intersections, rounded to the metre.
%   Street class: 1 = medium-density residential, 2 = high-density
%   residential, 3 = commercial (market) street.
W.name = {'A1','A2','A3','A4','B1','B2','B3','B4', ...
          'C1','C2','C3','C4','D1','D2','D3','D4'};
W.xy = [   0    0;  165   -8;  335    5;  505   -5;   % row A (south)
          -6  152;  158  160;  328  146;  497  155;   % row B (market street)
           4  305;  172  312;  340  298;  512  306;   % row C
          -3  452;  163  461;  331  447;  508  455];  % row D (north)
E = [ 1  2;  2  3;  3  4;      % row A
      5  6;  6  7;  7  8;      % row B
      9 10; 10 11; 11 12;      % row C
     13 14; 14 15; 15 16;      % row D
      1  5;  5  9;  9 13;      % column 1
      2  6;  6 10; 10 14;      % column 2
      3  7;  7 11; 11 15;      % column 3
      4  8;  8 12; 12 16;      % column 4
      2  5;                    % cut-through lane A2-B1
      6 11;                    % diagonal B2-C3
     12 15];                   % diagonal C4-D3
cls = [1 1 1, 3 3 3, 2 2 2, 1 1 1, 1 1 1, 2 2 2, 2 2 2, 1 1 1, 1 3 2]';
E = sort(E, 2);
W.E = E;
W.class = cls;
W.len = round(sqrt(sum((W.xy(E(:,1),:) - W.xy(E(:,2),:)).^2, 2)));
W.className = {'Residential (medium)', 'Residential (high)', 'Commercial'};
W.rho = [0.55 0.85 1.30];          % kg of kerbside waste per metre of street
W.demand = W.len .* W.rho(cls)';   % kg per street per collection day
% Transfer station S (truck base and dump), reached through two gates
W.gate = [1 4];                    % A1 and A4
W.gateLen = [1350 1200];           % access-road length to S (m)
end
