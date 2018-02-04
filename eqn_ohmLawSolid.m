function [g, Jrow] = eqn_ohmLawSolid(j,k,P,x,cur)

il = x(P.idx_il:P.nx:P.nx*P.nj, :);
phis = x(P.idx_phis:P.nx:P.nx*P.nj, :);

% solid conductivity and its derivative wrt cl
[sigma, ~] = fun_sigma(P);

% init jacobians rows
a = zeros(1,P.nx);
b = zeros(1,P.nx);
d = zeros(1,P.nx);

if j==1
    hx = P.L_neg/P.n_neg;
    
    g = -il(j,k);
    a = [];
    b(P.idx_il) = 1;

elseif j < P.bnd_sep_neg % negative electrode
    hx = P.L_neg/P.n_neg;
    
    g = -(phis(j,k)-phis(j-1,k)) - hx*cur/sigma(j) + hx*(il(j,k)+il(j-1,k))/2/sigma(j);
    
    b(P.idx_phis) = 1;
    b(P.idx_il) = -1/2/sigma(j)*hx;
    
    a(P.idx_phis) = -1;
    a(P.idx_il) = -1/2/sigma(j)*hx;
    
elseif j == P.bnd_sep_neg % set constant current at separtor il=cur
    hx = P.L_neg/P.n_neg;
    
    g = -(phis(j,k)-phis(j-1,k)) - hx*cur/sigma(j) + hx*(il(j,k)+il(j-1,k))/2/sigma(j);
    
    b(P.idx_phis) = 1;
    b(P.idx_il) = -1/2/sigma(j)*hx;
    
    a(P.idx_phis) = -1;
    a(P.idx_il) = -1/2/sigma(j)*hx;
    
elseif j < P.bnd_pos_sep % separator set phis = 0
    g = -phis(j,k);
    b(P.idx_phis) = 1;
    
elseif j < P.nj % positive electrode
    hx = P.L_pos/P.n_pos;
    
    g = -(phis(j+1,k)-phis(j,k)) - hx*cur/sigma(j) + hx*(il(j,k)+il(j+1,k))/2/sigma(j);
    
    b(P.idx_phis) = -1;
    b(P.idx_il) = -1/sigma(j)*hx/2;
    
    d(P.idx_phis) = 1;
    d(P.idx_il) = -1/sigma(j)*hx/2;
    
else %P.nj boudary il = 0
    g = -il(j,k);
    
    b(P.idx_il) = 1;
    d = [];
end

% concatenation of the the row vectors for the band jacobian matrix
Jrow = [a b d];