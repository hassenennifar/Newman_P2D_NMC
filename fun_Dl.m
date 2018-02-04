function [Dl, dDl] = fun_Dl(cl, P)
% cl in mol/m3
T = P.T; % in K

Dl = 10.^(-0.22*cl/1e3 - 4.43 - 54./(T-229-5*cl/1e3));

dDl = -10.^(54./(cl/200 - T + 229) - (11*cl)/50000 ...
    - 443/100)*log(10).*(27./(100*(cl/200 - T + 229).^2) + 11/50000);

Dl = Dl*0.01;
dDl = dDl*0.01;