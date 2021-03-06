function  dudt = dEeqdT_neg(sto)

dudt = 1e-3.*(0.00527 + 3.29927.*sto - 91.79326.*sto.^2 ...
       + 1004.91101.*sto.^3 - 5812.27813.*sto.^4 + 19329.75490.*sto.^5 ...
       - 37147.89470.*sto.^6 + 38379.18127.*sto.^7 - 16515.05308.*sto.^8) ...
       ./(1 - 48.09287.*sto + 1017.23480.*sto.^2 - 10481.80419.*sto.^3 ...
       + 59431.30001.*sto.^4 - 195881.64880.*sto.^5 + 374577.31520.*sto.^6 ...
       - 385821.16070.*sto.^7 + 165705.85970.*sto.^8);