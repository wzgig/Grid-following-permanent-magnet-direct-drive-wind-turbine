clear;
clc;

run("tuoMWRef_and_Para.m");
eqn1 = (Te - T_m - D_m*omega_m) / J == 0;
eqn2 = u_sq - R_s*i_sq - omega_e*Psi_sd == 0;
eqn3 = u_sd - R_s*i_sd + omega_e*Psi_sq == 0;
% eqn2 = (u_sq - R_s*i_sq - omega_e*L_sd*i_sd - omega_e*Psi_f)/L_sq == 0;
% eqn3 = (u_sd - R_s*i_sd + omega_e*L_sq*i_sq)/L_sd == 0;
eqn4 = i_sdref - i_sd == 0;
eqn5 = i_sqref - i_sq == 0;
eqn6 = n_ref - n == 0;
eqn7 = (u_sdref - u_sd)/T_d == 0;
eqn8 = (u_sqref - u_sq)/T_d == 0;
eqn9 = v_g_q == 0;
eqn10 = Kp_PLL*v_g_q+Ki_PLL*PLL_int+w_g == 0;
eqn11 = (P_s - P_dc)/(C_dc*U_dc) == 0;
eqn12 = U_dc - U_dcref == 0;
eqn13 = (v_g_d - u_gd - R_g*i_gd + w_g*L_g*i_gq)/L_g == 0;
eqn14 = (v_g_q - u_gq - R_g*i_gq - w_g*L_g*i_gd)/L_g == 0;
eqn15 = i_gdref - i_gd == 0;
eqn16 = i_gqref - i_gq == 0;
eqn17 = i_gD - iD == 0;
eqn18 = i_gQ - iQ == 0;

% disp(eqn1);disp(eqn2);disp(eqn3);disp(eqn4);disp(eqn5);disp(eqn6);disp(eqn7);disp(eqn8);
% disp(eqn9);disp(eqn10);disp(eqn11);disp(eqn12);disp(eqn13);disp(eqn14);disp(eqn15);disp(eqn16);
% disp(eqn17);disp(eqn18);

% syms omega_m Psi_sq Psi_sd Id_stator_int Iq_stator_int Speed_int u_sd u_sq ...
%     PLL_int thetapll U_dc Udc_int i_gd i_gq Id_grid_int Iq_grid_int;
[omega_m, Psi_sq, Psi_sd, Id_stator_int, Iq_stator_int, Speed_int,...
    u_sd, u_sq, PLL_int, thetapll, U_dc, Udc_int, i_gd, i_gq, Id_grid_int, Iq_grid_int,...
    v_w, i_gqref]...
                      = vpasolve(eval([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6, ...
                      eqn7,eqn8,eqn9,eqn10,eqn11,eqn12,eqn13,eqn14,eqn15,eqn16, ...
                      eqn17,eqn18]),...
                      [omega_m Psi_sq Psi_sd Id_stator_int Iq_stator_int Speed_int ...
                      u_sd u_sq PLL_int thetapll U_dc Udc_int i_gd i_gq Id_grid_int Iq_grid_int ...
                      v_w i_gqref],...
                      [-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf; ...
                      -inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf;-inf,inf; ...
                      -inf,inf;-inf,inf]);

inistate = [omega_m Psi_sq Psi_sd Id_stator_int Iq_stator_int Speed_int...
    u_sd u_sq PLL_int thetapll U_dc Udc_int i_gd i_gq Id_grid_int Iq_grid_int]'
state_size = size(inistate,1);

v_w = double(v_w);
i_gqref = double(i_gqref);
disp(v_w);disp(i_gqref);

dt = 0.00001;
tspan=dt:dt:10;
options = odeset('RelTol',1e-12,'AbsTol',...
    1e-12*ones(1,state_size));
[t,x] = ode45(@(t,x) tuoMW(t,x,v_w, i_gqref),tspan,double(inistate),options);

% run("tuoMWRef_and_Para.m");

%%%%%%%%%% states

omega_m = x(:,1);
Psi_sq = x(:,2);
Psi_sd = x(:,3);
Id_stator_int = x(:,4);  % 原x1：定子d轴电流环积分状态
Iq_stator_int = x(:,5);  % 原x2：定子q轴电流环积分状态
Speed_int = x(:,6);      % 原y1：转速环积分状态
u_sd = x(:,7);
u_sq = x(:,8);
PLL_int = x(:,9);        % 原xpll：锁相环积分状态
thetapll = x(:,10);
U_dc = x(:,11);
Udc_int = x(:,12);       % 原z1：直流电压环积分状态
i_gd = x(:,13);
i_gq = x(:,14);
Id_grid_int = x(:,15);   % 原z2：网侧d轴电流环积分状态
Iq_grid_int = x(:,16);   % 原z3：网侧q轴电流环积分状态

% plot(tspan, eval(i_sd), 'DisplayName', 'i_sd');
% hold on;
% plot(tspan, eval(i_sq), 'DisplayName', 'i_sq');
% hold on;
% plot(tspan, u_sd, 'DisplayName', 'u_sd');
% hold on;
% plot(tspan, u_sq, 'DisplayName', 'u_sq');
% hold on;
% plot(tspan, eval(Cp), 'DisplayName', 'Cp');
% hold on;
% plot(tspan, eval(T_m), 'DisplayName', 'T_m');
% hold on;
% % plot(tspan, eval(Tem), 'DisplayName', 'Tem');
% % hold on;
% plot(tspan, eval(Te), 'DisplayName', 'Te');
% hold on;
% plot(tspan, eval(omega_e), 'DisplayName', 'omega_e');
% hold on;
% % plot(tspan, eval(Tem_cmd), 'DisplayName', 'Tem_cmd');
% % hold on;
% plot(tspan, omega_m, 'DisplayName', 'omega_m');
% hold on;
figure(1);
plot(tspan, eval(P_m)/Sbase, 'DisplayName', 'P_m');
hold on;
plot(tspan, eval(-P_s)/Sbase, 'DisplayName', 'P_s');
hold on;
% plot(tspan, eval(-P_e)/Sbase, 'DisplayName', 'P_e');
legend('Location', 'best', 'Interpreter', 'tex');
% plot(tspan, eval(Q_s), 'DisplayName', 'Q_s');
% hold on;
% plot(tspan, eval(P_e), 'DisplayName', 'P_e');
% hold on;
% plot(tspan, eval(n), 'DisplayName', 'n');
% hold on;
% plot(tspan, eval(P_g), 'DisplayName', 'P_g');
% hold on;
% plot(tspan, eval(Q_g), 'DisplayName', 'Q_g');
% hold on;
% plot(tspan, eval(P_dc), 'DisplayName', 'P_dc');
% hold on;
% plot(tspan, eval(u_gd), 'DisplayName', 'u_gd');
% hold on;
% plot(tspan, eval(u_gq), 'DisplayName', 'u_gq');
% hold on;
% plot(tspan, i_gd, 'DisplayName', 'i_gd');
% hold on;
% plot(tspan, i_gq, 'DisplayName', 'i_gq');
% hold on;

figure(2);
plot(tspan, eval(P_g/Sbase), 'DisplayName', 'P_g^{pu}');
hold on;
plot(tspan, eval(Q_g/Sbase), 'DisplayName', 'P_g^{pu}');
hold on;
legend('Location', 'best', 'Interpreter', 'tex');   % Location参数指定图例位置，'best'表示自动选择最佳位置