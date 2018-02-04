function tw = mass(P)

% mass of positive electrode
m_pos = P.L_pos*(P.re*P.epsl_pos+P.rs_pos*P.epss_pos+P.rf*P.epf_pos);
% mass of separator
m_sep = P.L_neg*(P.re*P.epsl_sep+P.rc*(1-P.epsl_sep));

% mass of negative electrode
m_neg = P.L_neg*(P.re*P.epsl_neg+P.rs_neg*P.epss_neg+P.rf*P.epf_neg);

% mass of current collectors
cc1 = P.rcn*P.hcn+P.rcp*P.hcp;

tw = m_pos+m_sep+m_neg+cc1;