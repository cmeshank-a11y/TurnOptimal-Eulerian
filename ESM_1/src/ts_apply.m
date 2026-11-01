function partner = ts_apply(partner, Q)
%TS_APPLY  Write the pairing Q (k-by-2 half-edge pairs) into partner.
partner(Q(:,1)) = Q(:,2);
partner(Q(:,2)) = Q(:,1);
end
