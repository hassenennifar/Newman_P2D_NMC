function [g, Jrow] = eqn_curDivLiquid(j,k,P,x,cur)

il = x(P.idx_il:P.nx:P.nx*P.nj, :);
jn = x(P.idx_jn:P.nx:P.nx*P.nj, :);
phil = x(P.idx_phil:P.nx:P.nx*P.nj, :);
phis = x(P.idx_phis:P.nx:P.nx*P.nj, :);

% init jacobians rows
a = zeros(1,P.nx);
b = zeros(1,P.nx);
d = zeros(1,P.nx);

if j == 1  % il = 0
    hx = P.L_neg/P.n_neg;
    area = 3*P.epss_neg/P.Rp_neg;
    Cdl = P.Cdl_neg;
%     area = area/2;

    g = (il(j+1,k)-il(j,k))/hx/area/P.F - (1*jn(j,k)+jn(j+1,k))/2;
    
    a = [];
    
    b(P.idx_il) = 1/hx/area/P.F;
    b(P.idx_jn) = 1/2;
    
    if P.dt > 0
        g = g - Cdl/P.F/P.dt*(phis(j,k)-phil(j,k)-(phis(j,k-1)-phil(j,k-1)));
        b(P.idx_phil) = -Cdl/P.F/P.dt;
        b(P.idx_phis) = Cdl/P.F/P.dt;
    end
    

    d(P.idx_il) = -1/hx/area/P.F;
    d(P.idx_jn) = 1/2;
%     g = -il(j,k);
%     b(P.idx_il) = 1;

elseif j < P.bnd_sep_neg % negative electrode
    hx = P.L_neg/P.n_neg;
    area = 3*P.epss_neg/P.Rp_neg;
    Cdl = P.Cdl_neg;
    if j == P.bnd_sep_neg
%         area = area/2;
    end
    
    g = (il(j+1,k)-il(j,k))/hx/area/P.F - (1*jn(j,k)+jn(j+1,k))/2;
    
    b(P.idx_il) = 1/hx/area/P.F;
    b(P.idx_jn) = 1/2;
    
    if P.dt > 0
        g = g - Cdl/P.F/P.dt*(phis(j,k)-phil(j,k)-(phis(j,k-1)-phil(j,k-1)));
        b(P.idx_phil) = -Cdl/P.F/P.dt;
        b(P.idx_phis) = Cdl/P.F/P.dt;
    end
    
    d(P.idx_il) = -1/hx/area/P.F;
    d(P.idx_jn) = 1/2;

elseif j == P.bnd_sep_neg
    g = cur - il(j,k);
    b(P.idx_il) = 1;
    
elseif j <= P.bnd_pos_sep % separator 
    g = - il(j,k) + cur;
    
    b(P.idx_il) = 1;
%     a(P.idx_il) = -1;

% elseif j == P.bnd_pos_sep 
%     g = - il(j,k) + cur;
% 
%     b(P.idx_il) = 1;
    
elseif j <= P.nj % positive electrode
    hx = P.L_pos/P.n_pos;
    area = 3*P.epss_pos/P.Rp_pos;
    Cdl = P.Cdl_pos;
    if j == P.bnd_pos_sep || j == P.nj
%         area = area/2;
        %il(j,k) = il(j,k)/2;
    end
    
    g = (il(j,k)-il(j-1,k))/hx/area/P.F - (1*jn(j,k)+jn(j-1,k))/2;
    
    a(P.idx_il) = 1/hx/area/P.F; 
    a(P.idx_jn) = 1/2;
    
    b(P.idx_il) = -1/hx/area/P.F;
    b(P.idx_jn) = 1/2;
    
    if P.dt > 0
        g = g - Cdl/P.F/P.dt*(phis(j,k)-phil(j,k)-(phis(j,k-1)-phil(j,k-1)));
        b(P.idx_phil) = -Cdl/P.F/P.dt;
        b(P.idx_phis) = Cdl/P.F/P.dt;
    end
    
    if j == P.nj
        d = [];
    end
 
end
% concatenation of the the row vectors for the band jacobian matrix
Jrow = [a b d];