function [g0, g1] = Eeq_pos(sto, P)
    
g0 = 6.0826 - 6.9922*sto + 7.1062*sto.^2 - 0.54549e-4*exp(124.23*sto-114.2593) ...
    - 2.5947*sto.^3;
   
g1 = (35531*sto)/2500 - (100005293359112539557*exp((12423*sto)/100 ...
    - 2010070862904741/17592186044416))/14757395258967641292800 ...
    - (77841*sto.^2)/10000 - 34961/5000;
g1 = g1/P.csmax_pos;