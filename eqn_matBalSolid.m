function [g, Jrow] = eqn_matBalSolid(j,k,P,x,ts)

% extraction of required states
jn = x(P.idx_jn:P.nx:P.nx*P.nj, :);
css = x(P.idx_css:P.nx:P.nx*P.nj, :);

% init jacobians rows
a = zeros(1,P.nx);
b = zeros(1,P.nx);
d = zeros(1,P.nx);

if j == 1;
    a = [];
end
if j == P.nj
    d = [];
end

if j <= P.bnd_sep_neg
    Rp = P.Rp_neg;
elseif j >= P.bnd_pos_sep
    Rp = P.Rp_pos;
else % separator not defined
    g = -css(j,k);
    b(P.idx_css) = 1;
    Jrow = [a b d];
    return;
end

if P.dt == 0
    g = -css(j,k) + css(j,k-1);
    
    b(P.idx_css) = 1;
    Jrow = [a b d];
    return;
end

% superposition integral method
% [a1, sumpa] = calca(j,kt,ts, x, P);
a1 = P.a1(j);
sumpa = P.sumpa(j);

% zero equation
g = -jn(j,k) - Rp*(sumpa + a1*(css(j,k)-css(j,k-1))/P.dt);

% jacobian row
b(P.idx_css) = Rp*a1/P.dt;
b(P.idx_jn) = 1;

% concatenation of the the row vectors for the band jacobian matrix
Jrow = [a b d];