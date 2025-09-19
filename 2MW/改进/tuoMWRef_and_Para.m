f    = 60;
Fnom    = f;
Np  = 48;
w_g = 2*pi*f;
Sbase = 2e6;
Pnom = Sbase;
Vbase = 690;
Ibase = Sbase/(sqrt(3)*Vbase);
wbase = w_g;
Zbase = Vbase^2/Sbase;
Lbase = Zbase/w_g;
Cbase = 1/(wbase*Zbase);
Tbase = Sbase/(w_g/Np);
U_dcref = 1150;

R_g = 0.05;
L_g = 7e-3;
C_dc = 1e-3;

% 发电机参数
Psi_f = 3.88889;
L_sd = 1.8e-3;
L_sq = 1.8e-3;
R_s = 0.0026;
beta = 0;
pitch = beta;
J = 35000;
D_m = 0.078;
R_t = 36.6;

% **解决方案3: 调整控制器参数**
% 原来的PLL参数过于激进，降低增益
Kp_PLL = 100;    % 原来是350，降低到100
Ki_PLL = 1000;   % 原来是3200，降低到1000

% 适当调整其他控制器参数以确保稳定性
Kp_Id_stator = 0.8;   % 原来是1，稍微降低
Ki_Id_stator = 10;    % 原来是12，稍微降低
Kp_Iq_stator = 0.8;   % 原来是1，稍微降低
Ki_Iq_stator = 10;    % 原来是12，稍微降低
Kp_Speed = 80;        % 原来是100，稍微降低
Ki_Speed = 180;       % 原来是220，稍微降低
Kp_Udc = 0.8;         % 原来是1，稍微降低
Ki_Udc = 8;           % 原来是10，稍微降低
Kp_Id_grid = 0.8;     % 原来是1，稍微降低
Ki_Id_grid = 12;      % 原来是15，稍微降低
Kp_Iq_grid = 0.8;     % 原来是1，稍微降低
Ki_Iq_grid = 12;      % 原来是15，稍微降低

% 适当增大时间常数以降低系统带宽
T_d = 1/4000;         % 原来是1/6000，适当增大
T_m = 1/2000;         % 原来是1/3000，适当增大
T_trq = 60;

% **解决方案4: 调整xi的选择策略**
% 对于敏感的xi范围，使用更保守的设置
% 可以根据需要调整这个值
xi_target = -pi/6;    % 目标相角

% 检查xi是否在敏感范围内
if xi_target >= -pi/2 && xi_target <= -pi/6
    fprintf('警告: xi在敏感范围内 [%.2f, %.2f]\n', -pi/2, -pi/6);
    fprintf('已应用保守的控制器参数设置\n');
    
    % 在敏感范围内进一步降低PLL增益
    Kp_PLL = Kp_PLL * 0.7;
    Ki_PLL = Ki_PLL * 0.7;
end

% 目标有功功率和无功功率
P = -1;     % 发电模式下为负值，有功功率（p.u.）
Q = 0;      % 无功功率（p.u.）
V = 1;      % 母线电压幅值（p.u.）
xi = xi_target;  % 使用目标相角

% 计算 αβ 坐标系下的电流参考值
S_D0 = P/V;
S_Q0 = Q/V;
S_ab0 = S_D0 + 1j*S_Q0;
i_abs = abs(S_ab0);
ui_argdiff = angle(S_ab0);
i_arg = xi - ui_argdiff;
i_ab = i_abs * exp(1i * i_arg);
iQ = imag(i_ab);
iD = real(i_ab);

v_ab = V * exp(1i * xi);
v_a = real(v_ab);
v_b = imag(v_ab);

% 显示关键参数
fprintf('关键参数设置:\n');
fprintf('xi = %.4f rad (%.1f°)\n', xi, xi*180/pi);
fprintf('v_a = %.4f, v_b = %.4f\n', v_a, v_b);
fprintf('iD = %.4f, iQ = %.4f\n', iD, iQ);
fprintf('PLL增益: Kp = %.1f, Ki = %.1f\n', Kp_PLL, Ki_PLL);

% 定义符号变量用于稳态计算
syms omega_m Psi_sq Psi_sd Id_stator_int Iq_stator_int Speed_int...
    u_sd u_sq PLL_int thetapll U_dc Udc_int i_gd i_gq Id_grid_int Iq_grid_int;
syms v_w i_gqref

% 计算参考值和中间变量
omega_best = 8.1*v_w/R_t;
n_ref = omega_best*30/pi;

n = omega_m * 30/pi;
omega_e = omega_m*Np;
omega_r = omega_m * (Np/(Fnom*2*pi));

i_sd = (Psi_sd - Psi_f)/L_sd;
i_sq = Psi_sq/L_sq;

Te = sqrt(3)*Np*(Psi_f*i_sq);
Tem = Te * (1/(Np*Pnom/(2*pi*Fnom)));

% 风机模型
u1 = 0.5*pi;
u2 = 1.225;
lambda = (R_t*omega_m)/(v_w);
lambda_i = 1/(1/(lambda+0.08*beta)-0.035/(beta^3+1));
Cp = 0.51763*(116/lambda_i-0.4*beta-5)*exp(-21/lambda_i)+0.006795*lambda;
P_m = u1*u2*(R_t^2)*(v_w^3)*Cp;
T_m = (-1)*(P_m)/omega_m;

% 控制器
i_sdref = 0;
i_sqref = Kp_Speed*(n_ref - n) + Ki_Speed*Speed_int;
u_sqref = Kp_Iq_stator*(i_sqref - i_sq)+Ki_Iq_stator*Iq_stator_int+omega_e*Psi_f+omega_e*L_sd*i_sd;
u_sdref = Kp_Id_stator*(i_sdref - i_sd)+Ki_Id_stator*Id_stator_int-omega_e*L_sq*i_sq;

P_s = sqrt(3)*(u_sd * i_sd + u_sq * i_sq);
Q_s = sqrt(3)*(u_sq * i_sd - u_sd * i_sq);
P_e = Te * omega_m;

% 坐标变换和网侧控制
theta_g = thetapll;
Tg_DQdq = [cos(theta_g) -sin(theta_g); sin(theta_g) cos(theta_g)];
vg_dq_real = Vbase * Tg_DQdq * [v_a; v_b];
v_g_d = vg_dq_real(1);
v_g_q = -vg_dq_real(2);

Tg_dqDQ = [cos(-theta_g) -sin(-theta_g); sin(-theta_g) cos(-theta_g)];
ig_dq = (1/Ibase)*Tg_dqDQ*[i_gd;i_gq];

i_gD = ig_dq(1);
i_gQ = ig_dq(2);

i_gdref = Kp_Udc*(U_dc - U_dcref) + Ki_Udc*Udc_int;

u_gd = - Kp_Id_grid*(i_gdref - i_gd) - Ki_Id_grid*Id_grid_int - R_g*i_gd + v_g_d+i_gq*w_g*L_g;
u_gq = - Kp_Iq_grid*(i_gqref - i_gq) - Ki_Iq_grid*Iq_grid_int - R_g*i_gq + v_g_q-i_gd*w_g*L_g;

P_dc = sqrt(3)*(u_gd* i_gd + u_gq * i_gq);
P_g = sqrt(3)*(v_g_d * i_gd + v_g_q * i_gq);
Q_g = sqrt(3)*(v_g_q * i_gd - v_g_d * i_gq);