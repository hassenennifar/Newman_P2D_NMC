function [g, Jrow] = eqn_ohmLawLiquid(j,k,P,x)

cl = x(P.idx_cl:P.nx:P.nx*P.nj, :);
il = x(P.idx_il:P.nx:P.nx*P.nj, :);
phil = x(P.idx_phil:P.nx:P.nx*P.nj, :);

% liquid conductivity and its derivative wrt cl
[kappa, dkappa] = fun_kappa(cl(:,k), P);

% tranport number and its derivative wrt cl
[tp, dtp] = fun_tp(cl(:,k), P);

% liquid activity coefficient and its derivative wrt cl
[dlnf, d2lnf] = fun_dlnf(cl(:,k), P);

RTF2 = 2*P.R*P.T/P.F;

% init jacobians rows
a = zeros(1,P.nx);
b = zeros(1,P.nx);
d = zeros(1,P.nx);

if j == P.nj % boundary cond for phil = 0 at pos current collector
    
    g = -phil(j,k);
    
    d = [];
    b(P.idx_phil) = 1;
    
else
    if j >= P.bnd_pos_sep % positive electrode
        hx = P.L_pos/P.n_pos;
        epsl = P.epsl_pos;
    elseif j >= P.bnd_sep_neg % separator
        hx = P.L_sep/P.n_sep;
        epsl = P.epsl_sep;
    else % negative electrode
        hx = P.L_neg/P.n_neg;
        epsl = P.epsl_neg;
    end
    
    if j == 1
        a = [];
    end
    
    g = (phil(j+1,k)-phil(j,k))/hx + 1/epsl^1.5*il(j,k)/kappa(j)/2+1/epsl^1.5*il(j+1,k)/kappa(j+1)/2 ...
        - RTF2*( (1/cl(j,k)+ dlnf(j))*(1-tp(j)) + (1/cl(j+1,k) + dlnf(j+1))*(1-tp(j+1)))*(cl(j+1,k)-cl(j,k))/hx/2;
    
    b(P.idx_cl) = -1/epsl^1.5*il(j,k)*dkappa(j)/kappa(j)^2/2 ...
        - RTF2*(-1/cl(j,k)^2 + d2lnf(j)/2) *(1-tp(j))*(cl(j+1,k)-cl(j,k))/hx/2 ...
        + RTF2*(1/cl(j,k)+ dlnf(j))* dtp(j)/2 * (cl(j+1,k)-cl(j,k))/hx/2 ...
        + RTF2*((1/cl(j,k)+ dlnf(j))*(1-tp(j)) + (1/cl(j+1,k) + dlnf(j+1))*(1-tp(j+1)))/hx/2;
    
    b(P.idx_il) = 1/epsl^1.5/kappa(j)/2;
    
    b(P.idx_phil) = -1/hx;

    d(P.idx_il) = 1/epsl^1.5/kappa(j+1)/2;
    d(P.idx_cl) = -1/epsl^1.5*il(j+1,k)*dkappa(j+1)/kappa(j+1)^2/2 ...
        - RTF2*(-1/cl(j+1,k)^2 + d2lnf(j+1)/2) *(1-tp(j+1))*(cl(j+1,k)-cl(j,k))/hx/2 ...
        + RTF2*(1/cl(j+1,k)+ dlnf(j+1))* dtp(j+1)/2 * (cl(j+1,k)-cl(j,k))/hx/2 ...
        - RTF2*((1/cl(j,k)+ dlnf(j))*(1-tp(j)) + (1/cl(j+1,k) + dlnf(j+1))*(1-tp(j+1)))/hx/2;
    
    d(P.idx_phil) = 1/hx;
    
    b=-b*hx;%/(P.R*P.T/P.F);
    d=-d*hx;%/(P.R*P.T/P.F);
    g = g*hx;%/(P.R*P.T/P.F);
end

% concatenation of the the row vectors for the band jacobian matrix
Jrow = [a b d];