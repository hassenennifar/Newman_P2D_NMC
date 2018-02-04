function P = init_param()

P.F = 96485;  % Faraday constant (As/mol)
P.R = 8.314;  % gas constant (J/mol-K)

P.idx_cl = 1; % variable index liquid concentration 
P.idx_phil = 2; % variable index liquid potential
P.idx_css = 3; % variable index solid surface concentration
P.idx_il = 4; % variable index liquid current density
P.idx_jn = 5; % variable index pore-wall molar flux
P.idx_phis = 6; % variable index solid phase potential

P.nx = 6; % number of states

P.lim = 500;  % limit on number of iterations
P.tolrel = 1e-4; % relative tolerance
P.tolabs = 1e-10; % abolute tolerance

P.n_neg = 1; % number of nodes in negative electrode
P.n_sep = 1; % number of nodes in separator
P.n_pos = 1; % number of nodes in positive electrode
P.nj = P.n_neg+P.n_sep+P.n_pos+1; % total number of nodes

P.bnd_sep_neg = P.n_neg+1; % boundary node between separator and anode
P.bnd_pos_sep = P.n_neg+P.n_sep+1; % boundary node between cathode and separator

P.L_neg = 80e-6;  % thickness of negative electrode (m)
P.L_sep = 25e-6;  % thickness of separator (m)
P.L_pos = 68e-6;  % thickness of positive electrode (m)

P.hcn = 10e-06;  % thickness of negative electrode current collector (m)
P.hcp = 10e-06;  % thickness of positive electrode current collector (m)

P.cl0 = 1000;  % initial electrolyte concentration (mol/m3)
P.sto_neg = 0.8;  % initial stochiometric parameter for negative particle
P.sto_pos = 0.4;  % initial stochiometric parameter for positive particle
P.dtmax = 1.02;  % maximum time step size (s)

P.Ds_neg = 5e-14;  % diffusion coefficient in negative solid (m2/s)
P.Ds_pos = 3e-14;  % diffusion coefficient in positive solid (m2/s)

P.Ds_neg0 = P.Ds_neg; % save value at 25degC for Arrhenius coefficient
P.Ds_pos0 = P.Ds_pos; % save value at 25degC for Arrhenius coefficient

P.Rp_neg = 11e-6;  % radius of negative particles (m)
P.Rp_pos =  5e-6;  % radius of positive particles (m)

P.epsl_neg = 0.3;  %  volume fraction of electrolyte in negative electrode
P.epsl_sep = 0.45;  %  volume fraction of electrolyte in separator
P.epsl_pos = 0.3;  %  volume fraction of electrolyte in positive electrode

P.epf_neg = 0.1;  %  volume fraction of inert filler in negative electrode
P.epf_pos = 0.1;  %  volume fraction of inert filler in positive electrode

P.epss_neg = 1-P.epsl_neg-P.epf_neg; % volume fraction of solid phase in negative electrode
P.epss_pos = 1-P.epsl_pos-P.epf_pos; % volume fraction of solid phase in negative electrode

P.sigma_neg = 100;  %  conductivity of solid negative matrix (S/m)
P.sigma_pos = 0.5;  %  conductivity of solid positive matrix (S/m)

P.rka_neg = 2.22e-11;  %  reaction rate constant for negative insertion reaction
P.rka_pos = 6e-12;  %  reaction rate constant for positive insertion reaction

P.rka_neg0 = P.rka_neg; % save value at 25degC for Arrhenius coefficient
P.rka_pos0 = P.rka_pos; % save value at 25degC for Arrhenius coefficient

P.Rf_neg = 0.35e-2;  %  anode film resistance (out of place) (ohm-m2) 
P.Rf_pos = 0;  %  cathode film resistance (out of place) (ohm-m2)

P.Rf_neg0 = P.Rf_neg; % save value at 25degC for Arrhenius coefficient
P.Rf_pos0 = P.Rf_pos; % save value at 25degC for Arrhenius coefficient

P.cot_neg = 372;  %  coulombic capacity of negative material (mAh/g)
P.cot_pos = 274;  %  coulombic capacity of positive material (mAh/g)

P.re = 1324;  %  density of electrolyte (kg/m3)
P.rs_neg = 1800;  %  density of negative insertion material (kg/m3)
P.rs_pos = 5010;  %  density of positive insertion material (kg/m3)
P.rf = 1800;  %  density of inert filler (kg/m3)
P.rpl = 1780;  %  density of polymer phase (kg/m3)
P.rc = 552;  %  density of separator material (kg/m3)
P.rcn = 8954;  %  density of negative current collector (kg/m3)
P.rcp = 2707;  %  density of positive current collector (kg/m3)

P.csmax_neg = 31360; % 3600*P.cot_neg*P.rs_neg/P.F; % max solid phase concentration negative electrode (mol/m3)
P.csmax_pos = 52501; % 3600*P.cot_pos*P.rs_pos/P.F; % max solid phase concentration positive electrode (mol/m3)
P.cs0_neg = P.sto_neg*P.csmax_neg; % init solid phase concentration negative electrode (mol/m3)
P.cs0_pos = P.sto_pos*P.csmax_pos; % init solid phase concentration positive electrode (mol/m3)

P.tort_neg = 5.25;
P.tort_sep = 2.5;
P.tort_pos = 3.5;

P.Cdl_neg = 0.1; % F/m2
P.Cdl_pos = 0.1; % F/m2

P.htc = 5;  % heat-transfer coefficient with external medium (W/m2K)
P.Cp = 500;  % heat capacity of cell (J/kgK)
P.Tam  = 298.15;  % ambient temperature (K)
P.T  = P.Tam;  % temperature (K)

P.RG = 8e-4*0; % total resistance in foils, leads, and contacts, ohm-m2
P.RGext = P.RG/4; % resistance outside cell

P.EbarDs_neg = 0;  % activation energy, anode solid diffusion
P.EbarDs_pos = 0;  % activiation energy, cathode solid diffusion
P.Ebarkap = 4000;  % activation energy electrolyte conductivity
P.EbarD  = 4000;  % activation energy electrolyte diffusion
P.Ebarka_neg = 4000;  % activation energy negative kinetics
P.Ebarka_pos = 4000;  % activation energy positive kinetics
P.Ebarks1a = 4000;  % activation energy O2 side rxn
P.Ebarks1c = 4000;  % activation energy O2 side rxn
P.Ebarks2a = 4000;  % activation energy H2 side rxn
P.Ebarks2c = 4000;  % activation energy H2 side rxn
P.Ebarks3a = 4000;  % activation energy shuttle side rxn
P.Ebarks3c = 4000;  % activation energy shuttle side rxn
P.Ebarr_neg = 4000;  % activation energy, film resistance anode
P.Ebarr_pos = 4000;  % activation energy, film resistance cathode