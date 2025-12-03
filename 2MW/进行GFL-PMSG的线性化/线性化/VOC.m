function [dotx, yy] = VOC(t,x,u)

%%%%%%%%%%%%%%%%%%%%%%%%%% Rename %%%%%%%%%%%%%%%%%%%%%%%%
PLL_int       = x(1);
thetapll      = x(2);
U_dc          = x(3);
Udc_int       = x(4);
i_gd          = x(5);
i_gq          = x(6);
Id_grid_int   = x(7);
Iq_grid_int   = x(8);

%%%%%%%%%%%%%%%%%%%%%%%%%%  Paras  %%%%%%%%%%%%%%%%%%%%%%%%
v_a = u(1);   % the global voltage should be input
v_b = u(2);
P_s = u(3);
i_gqref = u(4);

U_dcref = 1150;

f   = 60;
w_g = 2*pi*f;

Sbase = 2e6;
Pnom = Sbase;
Vbase = 690;
Ibase = Sbase/(sqrt(3)*Vbase);
wbase = w_g;
Zbase = Vbase^2/Sbase;
Lbase = Zbase/w_g;
U_dcref = 1150;

Kp_Udc = 1;
Ki_Udc = 10;
Kp_Id_grid = 1;
Ki_Id_grid = 15;
Kp_Iq_grid = 1;
Ki_Iq_grid = 15;

Kp_PLL = 350;
Ki_PLL = 3200;
%滤波部分 
R_g = 0.05;
L_g = 7e-3;
C_dc = 1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%% Algebra %%%%%%%%%%%%%%%%%%%%%%%%
%% ========== 坐标变换（GSC的电压观测） ==========
% PLL
theta_g = thetapll;         % 网侧PLL估计的同步角
% 构造变换矩阵，将d-q轴量转换为定轴坐标系（用于测量或控制）
Tg_DQdq = [cos(theta_g) -sin(theta_g); sin(theta_g) cos(theta_g)];
vg_dq_real = Vbase * Tg_DQdq * [v_a; v_b];  % GSC中的实际电压
% 变换后各轴上的电压值提取
v_g_d = vg_dq_real(1);          % 网侧d轴电压
v_g_q = -vg_dq_real(2);         % 网侧q轴电压

i_gdref = Kp_Udc*(U_dc - U_dcref) + Ki_Udc*Udc_int;
% i_gqref = 0;
u_gd = - Kp_Id_grid*(i_gdref - i_gd) - Ki_Id_grid*Id_grid_int - R_g*i_gd + v_g_d+i_gq*w_g*L_g;
u_gq = - Kp_Iq_grid*(i_gqref - i_gq) - Ki_Iq_grid*Iq_grid_int - R_g*i_gq + v_g_q-i_gd*w_g*L_g;

% P_s = Pnom;
P_dc = sqrt(3)*(u_gd* i_gd + u_gq * i_gq);%(2-17)
P_g = sqrt(3)*(v_g_d * i_gd + v_g_q * i_gq);%(2-17)
Q_g = sqrt(3)*(v_g_q * i_gd - v_g_d * i_gq);%(2-18)

%%%%%%%%%%%%%%%%%%%%%%% Differential %%%%%%%%%%%%%%%%%%%%%
Tg_dqDQ = [cos(-theta_g) -sin(-theta_g); sin(-theta_g) cos(-theta_g)];
ig_dq = (1/Ibase)*Tg_dqDQ*[i_gd;i_gq];

i_gD = ig_dq(1);
i_gQ = ig_dq(2);

yy = [i_gD; i_gQ];


    dotx = [
         v_g_q;
         Kp_PLL*v_g_q+Ki_PLL*PLL_int+w_g;
         (P_s - P_dc)/(C_dc*U_dc);
         U_dc - U_dcref;
         (v_g_d - u_gd - R_g*i_gd + w_g*L_g*i_gq)/L_g;
         (v_g_q - u_gq - R_g*i_gq - w_g*L_g*i_gd)/L_g;
         i_gdref - i_gd;
         i_gqref - i_gq;
    ];
end


