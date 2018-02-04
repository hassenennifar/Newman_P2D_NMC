function x = init_states(P,cur)

% initial states
cl = P.cl0*ones(P.nj,1);
phil = zeros(P.nj,1);
il = zeros(P.nj,1);
css = [P.cs0_neg*ones(P.n_neg+1,1); zeros(P.n_sep-1,1); P.cs0_pos*ones(P.n_pos+1,1)];
jn = zeros(P.nj,1);

Ua = Eeq_neg(P.cs0_neg/P.csmax_neg, P);
Uc = Eeq_pos(P.cs0_pos/P.csmax_pos, P);
phis = [Ua*ones(P.n_neg+1,1); zeros(P.n_sep-1,1); Uc*ones(P.n_pos+1,1)];


% % linear spatial distribution of current density
% kappa = fun_kappa(cl, P);
% phil(P.nj) = 0;
% area = 3*P.epss_pos/P.Rp_pos;
% for j = P.nj:-1:P.bnd_pos_sep % positive
%     il(j) = cur*(P.nj-j)/P.n_pos;
%     jn(j) = -cur/P.F/P.L_pos/area;
%     if(j ~= P.nj)
%         phil(j) = phil(j+1)+P.L_pos/P.n_pos*il(j)/kappa(j); % ohmic drop is not accurate
%     end
%     phis(j) = Uc+phil(j)+jn(j)*P.F*P.Rf_pos;
% end
% 
% for j = P.bnd_pos_sep-1:-1:P.bnd_sep_neg+1 % separator
%     il(j) = cur;
%     jn(j) = 0;
%     phil(j) = phil(j+1)+P.L_sep/P.n_sep*il(j)/kappa(j); % ohmic drop is not accurate
%     phis(j) = 0;
% end
% 
% area = 3*P.epss_neg/P.Rp_neg;
% for j = P.bnd_sep_neg:-1:1 % negative electrode
%     il(j) = cur*(j-1)/P.n_neg;
%     jn(j) = cur/P.F/P.L_neg/area;
%     phil(j) = phil(j+1)+P.L_neg/P.n_neg*il(j)/kappa(j); % ohmic drop is not accurate
%     phis(j) = Ua*0+phil(j)+jn(j)*P.F*P.Rf_neg;
% end

% intial state vector
x = zeros(P.nx*P.nj, 1);
k = 1;
for i = 1:P.nx:P.nx*P.nj
    x(i+P.idx_cl-1, k) = cl((i-1)/P.nx+1);
    x(i+P.idx_phil-1, k) = phil((i-1)/P.nx+1);
    x(i+P.idx_css-1, k) = css((i-1)/P.nx+1);
    x(i+P.idx_il-1, k) = il((i-1)/P.nx+1);
    x(i+P.idx_jn-1, k) = jn((i-1)/P.nx+1);
    x(i+P.idx_phis-1, k) = phis((i-1)/P.nx+1);
end