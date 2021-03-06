function [mag, ang, imp] = run_eis(potenz, vinit, xinit, P)

k = 1;
ts(k) = 0;
I(k) = 0;
v(k) = vinit;
x(:,k) = xinit;

freq = 10^potenz; % in Hz
w = 2*pi*freq;
mesh = 20;
P.dt = 1/freq/mesh;

while ts(k) <= 40/freq
    k = k+1;
    ts(k) = ts(k-1)+P.dt;

    x(:,k) = x(:,k-1);

    curamp = .001;
    cur = curamp*sin(w*ts(k));
    I(k) = cur;

    [x, P, ~] = runModel(k,P,x,cur,ts);

    v(k) = terminalVoltage(k,x,P,cur);

%     utz(:,k) = zeros(P.nj,1);
%     [utz(:,k-1:k), qq(k), P] = temperature(x(:,k-1:k), utz(:,k-1:k), 2, cur, v(k), P);
%     T(k) = P.T-273.15;
end
phis = x(P.idx_phis:P.nx:P.nx*P.nj,:);
veis = phis(end,:)-phis(1,:);
fitobj = fit(ts(end-4*mesh+1:end)', ...
    veis(end-4*mesh+1:end)'-mean(veis(end-4*mesh+1:end)), 'sin1');

fitobj2 = fit(ts(end-4*mesh+1:end)', ...
    I(end-4*mesh+1:end)', 'sin1');


mag = fitobj.a1 / curamp;
ang = pi-fitobj.c1;
imp = mag * exp(1j*ang);

figure(10)
cla
plot(fitobj,ts, veis-mean(veis(end-4*mesh+1:end)))
