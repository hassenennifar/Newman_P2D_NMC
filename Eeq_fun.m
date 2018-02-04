function Eeq = Eeq_fun(x,P)
css = x(P.idx_css:P.nx:P.nx*P.nj);

Eeq = zeros(P.nj,1);
Eeq(1:P.bnd_sep_neg) = Eeq_neg(css(1:P.bnd_sep_neg)/P.csmax_neg,P);
Eeq(P.bnd_pos_sep:P.nj) = Eeq_pos(css(P.bnd_pos_sep:P.nj)/P.csmax_pos,P);