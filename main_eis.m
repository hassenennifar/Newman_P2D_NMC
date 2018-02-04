clear

vcut = 3; % stop condition for voltage
cur = 21.5*8; % applied current (A/m2)

%% Initialization

P = init_param;
P.dt = 0; % time step

x = init_states(P,cur);
xinit = x;

k = 1;
ts(k) = 0; % time 

v(k) = terminalVoltage(k,x,P,cur); % terminal voltage
qq(k) = 0;
T(k) = P.T-273.15;
utz = zeros(P.nj,k);
utz(1:P.bnd_sep_neg,k) = P.cs0_neg/P.csmax_neg;
utz(P.bnd_pos_sep:P.nj,k) = P.cs0_pos/P.csmax_pos;
% total lithium in salt
P.totLiold = 0.07;
I(k) = 0;

vinit = v(1);
% Icur = [zeros(100,1); 21.5*ones(1000, 1); zeros(500,1)]; % uncomment for dynamic current

%% run simulation

potenz = -4:.2:5;

for count = 1:length(potenz)
    disp(potenz(count))
    [mag, ang, imp] = run_eis(potenz(count), vinit, xinit, P);
    
    Zmag(count) = mag;
    Zang(count) = ang;
    Zimp(count) = imp;
    
%     
    figure(20)
    hold on
    cla
    plot(-Zmag .* exp(1j*(Zang+pi/20)), '^-b')
%     
%     pause(.01)

end