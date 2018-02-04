clear

vcut = 3; % stop condition for voltage
cur = 30*4; % applied current (A/m2)

%% Initialization

P = init_param;
P.dt = 0; % time step

x = init_states(P,cur);

k = 1;
ts(k) = 0; % time 

v(k) = terminalVoltage(k,x,P,cur); % terminal voltage
qq(k) = 0;
T(k) = P.T-273.15;
I(k) = cur;
utz = zeros(P.nj,k);
utz(1:P.bnd_sep_neg,k) = P.cs0_neg/P.csmax_neg;
utz(P.bnd_pos_sep:P.nj,k) = P.cs0_pos/P.csmax_pos;
Eeq(:,k) = Eeq_fun(x(:,k),P);
phis(:,k) = x(P.idx_phis:P.nx:P.nx*P.nj,k);
phil(:,k) = x(P.idx_phil:P.nx:P.nx*P.nj,k);
% total lithium in salt
P.totLiold = 0.07;

% Icur = [zeros(100,1); 21.5*ones(1000, 1); zeros(500,1)]; % uncomment for dynamic current

%% run simulation

stop_simulation = false;
stop_condition = 'v(k) < vcut';
stop_condition_reached = false;
plating_condition = false;
while ~stop_simulation
    k = k+1;
    ts(k) = ts(k-1)+P.dt;
    I(k) = cur;
    x(:,k) = x(:,k-1);
    
    tic
    [x, P, iter(k-1)] = runModel(k,P,x,cur,ts);
    titer(k-1) = toc;
    
    v(k) = terminalVoltage(k,x,P,cur);
    
    utz(:,k) = zeros(P.nj,1);
    [utz(:,k-1:k), qq(k), P] = temperature(x(:,k-1:k), utz(:,k-1:k), 2, cur, v(k), P);
    T(k) = P.T-273.15;
    Eeq(:,k) = Eeq_fun(x(:,k),P);
    
    phis(:,k) = x(P.idx_phis:P.nx:P.nx*P.nj,k);
    phil(:,k) = x(P.idx_phil:P.nx:P.nx*P.nj,k);
    if phis(P.bnd_sep_neg,k)-phil(P.bnd_sep_neg,k) < 0.001
        plating_condition = true;
        kplating = k;
        potold = (phis(P.bnd_sep_neg,1)-phil(P.bnd_sep_neg,1)-Eeq(P.bnd_sep_neg,1));
        curold = 0;
        pot = phis(P.bnd_sep_neg,k)-phil(P.bnd_sep_neg,k)-Eeq(P.bnd_sep_neg,k);
        Rz = (pot-potold)/(cur-curold);
        cur = (-Eeq(P.bnd_sep_neg,k)-potold)/Rz;
%         cur = cur/2;
    end
    

    if P.dt <= 0
        P.dt = 0.025;
    elseif ts(k) < 1
        P.dt = 0.025;
%     elseif plating_condition
%         P.dt = 0.0;
%         plating_condition = false;
    elseif ~stop_condition_reached
        P.dt = 5;
    end
    
    if eval(stop_condition)
        stop_condition_reached = true;
        k = k-1;
        P.dt = P.dt/2;
        if P.dt < 0.01;
            stop_simulation = true;
            k = k+1;
        end
    end
    
end

%%



figure
hold on
plot(ts, Eeq([ P.bnd_sep_neg],:), '.-')
plot(ts, phis([P.bnd_sep_neg],:), '.-')
plot(ts,phis(P.bnd_sep_neg,:)-phil(P.bnd_sep_neg,:), '.-')

figure
plot(ts, phis([1 P.bnd_sep_neg],:)-Eeq([1 P.bnd_sep_neg],:), '.')
%%
% figure
% pp = 0;
% LL_neg = [0; cumsum(repmat(P.L_neg/P.n_neg,P.n_neg,1))];
% LL_sep = [LL_neg(end); cumsum(repmat(P.L_sep/P.n_sep,P.n_sep,1))];
% LL_pos = [LL_sep(end); cumsum(repmat(P.L_pos/P.n_pos,P.n_pos,1))];
% timepoints = [2 100 300 k];
% for tt = timepoints
%     pp = pp+1;
%     subplot(length(timepoints),1,pp)
%     plot(LL_neg, phis(1:P.bnd_sep_neg,tt) -phil(1:P.bnd_sep_neg,tt), 'o-')
%     legend(['time = ' num2str(ts(tt)) ' s'])
% end
%%
figure, plot(ts,I)
figure, plot(ts,phis(P.bnd_sep_neg,:)-phil(P.bnd_sep_neg,:), '.')
