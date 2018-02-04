function dudt = dEeqdT_pos(sto)

dudt = (-0.19952+0.92837.*sto-1.36455.*sto.^2+0.61154.*sto.^3)./ ...
    (1-5.66148.*sto+11.47636.*sto.^2-9.82431.*sto.^3+3.04876.*sto.^4);
dudt = dudt./1000; %V./K