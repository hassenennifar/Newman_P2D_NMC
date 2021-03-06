function [g, Jrow] = eqn_matBalLiquid(j,k,P,x)

cl = x(P.idx_cl:P.nx:P.nx*P.nj, :);
il = x(P.idx_il:P.nx:P.nx*P.nj, :);
phil = x(P.idx_phil:P.nx:P.nx*P.nj, :);
phis = x(P.idx_phis:P.nx:P.nx*P.nj, :);

% diffusion coefficient and its derivative wrt cl
[Dl, dDl] = fun_Dl(cl(:,k), P);

% tranport number and its derivative wrt cl
[tp, dtp] = fun_tp(cl(:,k), P);

% init jacobians rows
a = zeros(1,P.nx);
b = zeros(1,P.nx);
d = zeros(1,P.nx);

if P.dt == 0
    if j == 1;
        a = [];
    end
    if j == P.nj
        d = [];
    end
    
    g = cl(j,k-1) - cl(j,k);
    
    b(P.idx_cl) = 1;
    Jrow = [a b d];
    return;
end

if j == 1 % left boundary at negative current collector
    hx = P.L_neg/P.n_neg;
    epsl = P.epsl_neg;
    area = 3*P.epss_neg/P.Rp_neg;
    Cdl = P.Cdl_neg;
    
    g = epsl*hx/P.dt/8*(cl(j+1,k)-cl(j+1,k-1)+3*cl(j,k)-3*cl(j,k-1)) ...
        - epsl^1.5*(Dl(j+1)+Dl(j))/2*(cl(j+1,k)-cl(j,k))/hx ...
        - (1-(tp(j)+tp(j+1))/2)*(il(j,k)+il(j+1,k))/P.F/2 ...
        - area*Cdl*tp(j)*hx/P.F/P.dt*(phis(j,k)-phil(j,k)-(phis(j,k-1)-phil(j,k-1)));
    
    % jacobian row at left node j-1 doesn't exist
    a = [];
    
    b(P.idx_cl) = epsl*hx/P.dt*3/8 + epsl^1.5*(Dl(j+1)+Dl(j))/2/hx ...
        - epsl^1.5*dDl(j)/2*(cl(j+1,k)-cl(j,k))/hx ...
        + dtp(j)/2*(il(j,k)+il(j+1,k))/P.F/2;
    
    b(P.idx_il) = -(1-(tp(j)+tp(j+1))/2)/P.F/2;
    b(P.idx_phil) = area*Cdl*tp(j)*hx/P.F/P.dt;
    b(P.idx_phis) = -area*Cdl*tp(j)*hx/P.F/P.dt;
    
    d(P.idx_cl) = epsl*hx/P.dt/8 - epsl^1.5*(Dl(j+1)+Dl(j))/2/hx ...
        - epsl^1.5*dDl(j+1)/2*(cl(j+1,k)-cl(j,k))/hx ...
        + dtp(j+1)/2*(il(j,k)+il(j+1,k))/P.F/2;
    d(P.idx_il) = -(1-(tp(j)+tp(j+1))/2)/P.F/2;
    
elseif j == P.nj % right boundary at positive current collector
    epsl = P.epsl_pos;
    hx = P.L_pos/P.n_pos;
    area = 3*P.epss_pos/P.Rp_pos;
    Cdl = P.Cdl_pos;
    
    g = epsl*hx/P.dt/8*(cl(j-1,k)-cl(j-1,k-1)+3*cl(j,k)-3*cl(j,k-1)) ...
        + epsl^1.5*(Dl(j)+Dl(j-1))/2*(cl(j,k)-cl(j-1,k))/hx ...
        + (1-(tp(j)+tp(j-1))/2)*(il(j,k)+il(j-1,k))/P.F/2 ...
        + area*Cdl*tp(j)*hx/P.F/P.dt*(phis(j,k)-phil(j,k)-(phis(j,k-1)-phil(j,k-1)));
    
    a(P.idx_cl) = epsl*hx/P.dt/8 - epsl^1.5*(Dl(j)+Dl(j-1))/2/hx ...
        + epsl^1.5*dDl(j-1)/2*(cl(j,k)-cl(j-1,k))/hx ...
        - dtp(j-1)/2*(il(j,k)+il(j-1,k))/P.F/2;
    
    a(P.idx_il) = (1-(tp(j)+tp(j-1))/2)/P.F/2;
    
    b(P.idx_cl) = epsl*hx/P.dt*3/8 + epsl^1.5*(Dl(j)+Dl(j-1))/2/hx ...
        + epsl^1.5*dDl(j)/2*(cl(j,k)-cl(j-1,k))/hx ...
        - dtp(j)/2*(il(j,k)+il(j-1,k))/P.F/2;
    b(P.idx_il) = (1-(tp(j)+tp(j-1))/2)/P.F/2;
    b(P.idx_phil) = -area*Cdl*tp(j)*hx/P.F/P.dt;
    b(P.idx_phis) = area*Cdl*tp(j)*hx/P.F/P.dt;
    
    d = [];

else
    
    if j <= P.bnd_sep_neg % negative electrode
        hx = P.L_neg/P.n_neg;
        epsl = P.epsl_neg;
        area = 3*P.epss_neg/P.Rp_neg;
        Cdl = P.Cdl_neg;
    elseif j > P.bnd_pos_sep % positive electrode
        hx = P.L_pos/P.n_pos;
        epsl = P.epsl_pos;
        area = 3*P.epss_pos/P.Rp_pos;
        Cdl = P.Cdl_pos;
    else
        hx = P.L_sep/P.n_sep;
        epsl = P.epsl_sep;
        area = 0;
        Cdl = 0;
    end
    
    % left side
    g = epsl*hx/P.dt/8*(cl(j-1,k)-cl(j-1,k-1)+3*cl(j,k)-3*cl(j,k-1)) ...
        + epsl^1.5*(Dl(j)+Dl(j-1))/2*(cl(j,k)-cl(j-1,k))/hx ...
        + (1-(tp(j)+tp(j-1))/2)*(il(j,k)+il(j-1,k))/P.F/2 ...
        + area*Cdl*tp(j)*hx/P.F/P.dt*(phis(j,k)-phil(j,k)-(phis(j,k-1)-phil(j,k-1)));
    
    a(P.idx_cl) = epsl*hx/P.dt/8 - epsl^1.5*(Dl(j)+Dl(j-1))/2/hx ...
        + epsl^1.5*dDl(j-1)/2*(cl(j,k)-cl(j-1,k))/hx ...
        - dtp(j-1)/2*(il(j,k)+il(j-1,k))/P.F/2;
    
    a(P.idx_il) = (1-(tp(j)+tp(j-1))/2)/P.F;
    
    b(P.idx_cl) = epsl*hx/P.dt*3/8 + epsl^1.5*(Dl(j)+Dl(j-1))/2/hx ...
        + epsl^1.5*dDl(j)/2*(cl(j,k)-cl(j-1,k))/hx ...
        - dtp(j)/2*(il(j,k)+il(j-1,k))/P.F/2;
    b(P.idx_il) = (1-(tp(j)+tp(j-1))/2)/P.F;
    b(P.idx_phil) = -area*Cdl*tp(j)*hx/P.F/P.dt;
    b(P.idx_phis) = area*Cdl*tp(j)*hx/P.F/P.dt;
    
    % right side
    if j < P.bnd_sep_neg % negative electrode
        hx = P.L_neg/P.n_neg;
        epsl = P.epsl_neg;
    elseif j >= P.bnd_pos_sep % positive electrode
        hx = P.L_pos/P.n_pos;
        epsl = P.epsl_pos;
    else
        hx = P.L_sep/P.n_sep;
        epsl = P.epsl_sep;
    end
    
    g = g + epsl*hx/P.dt/8*(cl(j+1,k)-cl(j+1,k-1)+3*cl(j,k)-3*cl(j,k-1)) ...
        - epsl^1.5*(Dl(j+1)+Dl(j))/2*(cl(j+1,k)-cl(j,k))/hx ...
        - (1-(tp(j)+tp(j+1))/2)*(il(j,k)+il(j+1,k))/P.F/2;
    
    b(P.idx_cl) = b(P.idx_cl) + epsl*hx/P.dt*3/8 + epsl^1.5*(Dl(j+1)+Dl(j))/2/hx ...
        - epsl^1.5*dDl(j)/2*(cl(j+1,k)-cl(j,k))/hx ...
        + dtp(j)/2*(il(j,k)+il(j+1,k))/P.F/2;
    
    b(P.idx_il) = b(P.idx_il) - (1-(tp(j)+tp(j+1))/2)/P.F/2 ;
    
    d(P.idx_cl) = epsl*hx/P.dt/8 - epsl^1.5*(Dl(j+1)+Dl(j))/2/hx ...
        - epsl^1.5*dDl(j+1)/2*(cl(j+1,k)-cl(j,k))/hx ...
        + dtp(j+1)/2*(il(j,k)+il(j+1,k))/P.F/2;
    d(P.idx_il) = - (1-(tp(j)+tp(j+1))/2)/P.F/2 ;
    
end

% concatenation of the the row vectors for the band jacobian matrix
g = -g;
Jrow = [a b d];

if find(isnan(Jrow))
    error('NaN detected');
end