function [utz, qq, P] = temperature(x, utz, k, cur, v, P)

% calculate utilization
jn = x(P.idx_jn:P.nx:P.nx*P.nj, :);
css = x(P.idx_css:P.nx:P.nx*P.nj, :);

area_neg = 3*P.epss_neg/P.Rp_neg;
utz(1:P.bnd_sep_neg,k) = utz(1:P.bnd_sep_neg,k-1) - ...
    P.dt*area_neg/P.epss_neg/P.csmax_neg*jn(1:P.bnd_sep_neg,k);

area_pos = 3*P.epss_pos/P.Rp_pos;
utz(P.bnd_pos_sep:P.nj,k) = utz(P.bnd_pos_sep:P.nj,k-1) - ...
    P.dt*area_pos/P.epss_pos/P.csmax_pos*jn(P.bnd_pos_sep:P.nj,k);


%% calculate equilibrium potential and entropy in both electrodes
Eeq = zeros(P.nj,1);
dEeqdT = zeros(P.nj,1);
% utz(1:P.bnd_sep_neg,k) = css(1:P.bnd_sep_neg,k)/P.csmax_neg;
Eeq(1:P.bnd_sep_neg) = Eeq_neg(utz(1:P.bnd_sep_neg,k), P);
dEeqdT(1:P.bnd_sep_neg) = dEeqdT_neg(utz(1:P.bnd_sep_neg,k));

% utz(P.bnd_pos_sep:P.nj,k) = css(P.bnd_pos_sep:P.nj,k)/P.csmax_pos;
Eeq(P.bnd_pos_sep:P.nj) = Eeq_pos(utz(P.bnd_pos_sep:P.nj,k), P);
dEeqdT(P.bnd_pos_sep:P.nj) = dEeqdT_pos(utz(P.bnd_pos_sep:P.nj,k));


%% calculate heat generation in both electrodes
qq1_neg = jn(1:P.bnd_sep_neg-1,k).*(Eeq(1:P.bnd_sep_neg-1)-P.T*dEeqdT(1:P.bnd_sep_neg-1,1));
qq2_neg = jn(2:P.bnd_sep_neg,k).*(Eeq(2:P.bnd_sep_neg)-P.T*dEeqdT(2:P.bnd_sep_neg,1));
qq_neg = area_neg*P.F*P.L_neg/P.n_neg*sum(qq1_neg+qq2_neg)/2;


qq1_pos = jn(P.bnd_pos_sep:P.nj-1,k).*(Eeq(P.bnd_pos_sep:P.nj-1)-P.T*dEeqdT(P.bnd_pos_sep:P.nj-1,1));
qq2_pos = jn(P.bnd_pos_sep+1:P.nj,k).*(Eeq(P.bnd_pos_sep+1:P.nj)-P.T*dEeqdT(P.bnd_pos_sep+1:P.nj,1));
qq_pos = area_pos*P.F*P.L_pos/P.n_pos*sum(qq1_pos+qq2_pos)/2;

%% calculate total heat generation
qq = -cur*v-(qq_neg+qq_pos)-cur^2*P.RGext;

%% calculate new temperature
tw = mass(P);
P.T = 1/(1+P.dt*P.htc/P.Cp/tw)*(P.T+(P.dt/(P.Cp*tw))*(P.htc*P.Tam+qq));

%% update parameters with Arrhenius factor
% update solid diffusion coefficients with new temperature
P.Ds_neg = P.Ds_neg0*exp((P.EbarDs_neg)*(P.T-298.15)/(P.T*298.15));
P.Ds_pos = P.Ds_pos0*exp((P.EbarDs_pos)*(P.T-298.15)/(P.T*298.15));

% update reaction rate constants with new temperature
ti0_neg = exp((P.Ebarka_neg)*(P.T-298.15)/(P.T*298.15));
P.rka_neg = P.rka_neg0*ti0_neg;

ti0_pos = exp((P.Ebarka_pos)*(P.T-298.15)/(P.T*298.15));
P.rka_pos = P.rka_pos0*ti0_pos;

% update film resistances with new temperature
P.Rf_neg = P.Rf_neg0*exp((P.Ebarr_neg)*(298.15-P.T)/(P.T*298.15)); % anode film resistance
P.Rf_pos = P.Rf_pos0*exp((P.Ebarr_pos)*(298.15-P.T)/(P.T*298.15)); % cathode film resistance
