function v = terminalVoltage(k,x,P,cur)

phis = x(P.idx_phis:P.nx:P.nx*P.nj, :);

if k == 1
    v = phis(end,k)-phis(1,k);
else
    v = phis(end,k)-phis(1,k) - P.RG*cur;
end