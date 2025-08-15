
clear;
clc;


run("Ref_and_Para.m");

eqn1 = (p/J)*(T_e - T_shaft) == 0; % state1
eqn2 = (1/(2*H_WT))*(T_wt - T_shaft) == 0;% state20
eqn3 = Turbine_w - w_m/p;  % state21
eqn4 = u_sd - Rs*I_sd + w_g*Psi_qs == 0; % state2
eqn5 = u_sq - Rs*I_sq - w_g*Psi_ds == 0; % state3
eqn6 = u_rd - Rr*I_rd + w_r*Psi_qr == 0; % state4
eqn7 = u_rq - Rr*I_rq - w_r*Psi_dr == 0; % state5

[Turbine_w w_m Psi_ds Psi_qs Psi_dr Psi_qr shaft_w]...
                      = vpasolve(eval([eqn1, eqn2, eqn3, eqn4, eqn5, eqn6, eqn7]),...
                      [Turbine_w w_m Psi_ds Psi_qs Psi_dr Psi_qr shaft_w],...
                      [40,300;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf]);

inistate = [w_m Psi_ds Psi_qs Psi_dr Psi_qr Turbine_w shaft_w]';
state_size = size(inistate,1);

dt = 0.00001;
tspan=dt:dt:5;
options = odeset('RelTol',1e-12,'AbsTol',...
    1e-12*ones(1,state_size));
[t,x] = ode45(@(t,x) DFIG(t,x),tspan,double(inistate),options);

run("Ref_and_Para.m");

%%%%%%%%%% states

w_m    = x(:,1);
Psi_ds = x(:,2);
Psi_qs = x(:,3);
Psi_dr = x(:,4);
Psi_qr = x(:,5);

% Turbine_w = x(:,2);
% shaft_w= x(:,3);
plot(tspan,eval(T_e)/Tbase);
hold on
plot(tspan,eval(Ptot)/Sbase);
plot(tspan,eval((u_rq.*I_rq + u_rd.*I_rd)))
plot(tspan,eval((u_rq.*I_rq + u_rd.*I_rd)/(u_sq.*I_sq + u_sd.*I_sd)));
hold on
plot(tspan,eval(w_r)/w_g);
plot(tspan,w_m);
plot(tspan,eval(I_rq*u)/Ibase);
hold on
plot(tspan,eval(I_rd*u)/Ibase);
plot(tspan,eval(I_sd));
plot(tspan,eval(I_sq));
hold on
plot(tspan,eval(I_rd));
plot(tspan,eval(I_rq));
plot(tspan,eval(u_rd));

plot(tspan,eval(u_rq));
