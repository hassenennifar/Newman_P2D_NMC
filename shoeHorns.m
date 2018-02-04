function [xknew, nconcflag] = shoeHorns(xknew,x,k,P)

clnew = xknew(P.idx_cl:P.nx:P.nx*P.nj);
philnew = xknew(P.idx_phil:P.nx:P.nx*P.nj);
cssnew = xknew(P.idx_css:P.nx:P.nx*P.nj);
jnnew = xknew(P.idx_jn:P.nx:P.nx*P.nj);
ilnew = xknew(P.idx_il:P.nx:P.nx*P.nj);
phisnew = xknew(P.idx_phis:P.nx:P.nx*P.nj);

cl = x(P.idx_cl:P.nx:P.nx*P.nj);
phil = x(P.idx_phil:P.nx:P.nx*P.nj);
css = x(P.idx_css:P.nx:P.nx*P.nj);
% jn = x(P.idx_jn:P.nx:P.nx*P.nj, :);
% il = x(P.idx_il:P.nx:P.nx*P.nj, :);
phis = x(P.idx_phis:P.nx:P.nx*P.nj);

nerr = 0;


for j = 1:P.nj

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  shoe horns:

    if(clnew(j) < cl(j,k)/100)
        clnew(j) = cl(j,k)/100;
        %  Concentration of lithium in solution in the positive is kept from going to zero by the following.
        %  The minimum value can be adjusted depending to the specific scenario being tested.
    end
    if(clnew(j) > cl(j,k)*100)
        clnew(j) = cl(j,k)*100;
    end
    if (philnew(j) < (phil(j,k)-0.02))
        philnew(j) = phil(j,k)-0.02;
    end
    if (philnew(j) > (phil(j,k)+0.02))
        philnew(j) = phil(j,k)+0.02;
    end
    if (phisnew(j) < (phis(j,k)-0.02))
        phisnew(j) = phis(j,k)-0.02;
    end
    if (phisnew(j) > (phis(j,k)+0.02))
        phisnew(j) = phis(j,k)+0.02;
    end
    if (philnew(j) >  9.9)
        philnew(j) =  9.9;
    end
    if (philnew(j) < -9.9)
        philnew(j) = -9.9;
    end
    if (phisnew(j) >  9.9)
        phisnew(j) =  9.9;
    end
    if (phisnew(j) < -9.9)
        phisnew(j) = -9.9;
    end

%     for mpa = 1:P.npa
        if (j  >=  P.bnd_pos_sep)
            if(cssnew(j) < css(j,k)/100)
                nerr = nerr+1;
                cssnew(j) = css(j,k)/100; %  use cs min
            end
            if(cssnew(j) > css(j,k)*100)
                cssnew(j) = css(j,k)*100;
            end
            if(P.csmax_pos-cssnew(j) <= (P.csmax_pos-css(j,k))/100)
                nerr = nerr+1;
                cssnew(j) = P.csmax_pos-(P.csmax_pos-css(j,k))/100;
            end
            if(cssnew(j) >= P.csmax_pos)
                cssnew(j) = 0.999999*P.csmax_pos;
            end
        elseif (j <= P.bnd_sep_neg)

            if (cssnew(j) < css(j,k)/100)
                nerr = nerr+1;
                cssnew(j) = css(j,k)/100; %  use cs min
            end
            if (cssnew(j) > css(j,k)*100)
                cssnew(j) = css(j,k)*100;
            end
            if (P.csmax_neg-cssnew(j) <= (P.csmax_neg-css(j,k))/100)
                nerr = nerr+1;
                cssnew(j) = P.csmax_neg-(P.csmax_neg-css(j,k))/100;
            end
            if (cssnew(j) >= P.csmax_neg)
                cssnew(j) = 0.999999*P.csmax_neg;
            end
        end

        %  Not shoe-horns, but trips to force smaller time steps
        if (j >= P.bnd_pos_sep || j <= P.bnd_sep_neg)
            if(cssnew(j) < 1.0e-10)
                cssnew(j)  =  1e-10;
            elseif(cssnew(j) < 1.0e-3)
                nconcflag  =  1;
            end
        end
%     end %  mpa

    if(clnew(j) < 1.0e-10)
        clnew(j)  =  1e-10;
        %   c(P.kj,j)  =  0
    end
    if(clnew(j) < 1.0e-03)
        nconcflag  =  2;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%     for i = 1:n+1
%         xx(i,j) = c(i,j); % label_25
%     end
end

%  check for convergence
%  normalize concentrations
%  calculate total amount of lithium in solution
totLinew = 0;
for j = 1:P.nj
    fh = P.L_pos*P.epsl_pos/P.n_pos;
    if (j <= P.bnd_pos_sep)
        fh = P.L_sep*P.epsl_sep/P.n_sep;
    end
    if (j <= P.bnd_sep_neg)
        fh = P.L_neg*P.epsl_neg/P.n_neg;
    end
    if (j == 1)
        fh = 0.5*P.L_neg*P.epsl_neg/P.n_neg;
    end
    if (j == P.bnd_sep_neg)
        fh = 0.5*(P.L_neg*P.epsl_neg/P.n_neg+P.L_sep*P.epsl_sep/P.n_sep);
    end
    if(j == P.bnd_pos_sep)
        fh = 0.5*(P.L_pos*P.epsl_pos/P.n_pos+P.L_sep*P.epsl_sep/P.n_sep);
    end
    if (j == P.nj)
        fh = 0.5*P.L_pos*P.epsl_pos/P.n_pos;
    end
    totLinew = totLinew + fh*cl(j,k);
end
factor = P.totLiold/totLinew;
for j = 1:P.nj
    cl(j,k) = factor*cl(j,k);
end

for i = 1:P.nx:P.nx*P.nj
    xknew(i+P.idx_cl-1) = clnew((i-1)/P.nx+1);
    xknew(i+P.idx_phil-1) = philnew((i-1)/P.nx+1);
    xknew(i+P.idx_css-1) = cssnew((i-1)/P.nx+1);
    xknew(i+P.idx_il-1) = ilnew((i-1)/P.nx+1);
    xknew(i+P.idx_jn-1) = jnnew((i-1)/P.nx+1);
    xknew(i+P.idx_phis-1) = phisnew((i-1)/P.nx+1);
end