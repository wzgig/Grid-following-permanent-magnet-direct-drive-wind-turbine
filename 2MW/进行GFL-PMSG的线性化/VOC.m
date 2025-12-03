function [dotx, yy] = VOC(t,x,u)

%%%%%%%%%%%%%%%%%%%%%%%%%% Rename %%%%%%%%%%%%%%%%%%%%%%%%
i_gd = x(1);
i_gq = x(2);
v_ref = x(3);
theta = x(4);
Phi_vd = x(5);
Phi_vq = x(6);
Phi_cd = x(7);
Phi_cq = x(8);

%%%%%%%%%%%%%%%%%%%%%%%%%%  Paras  %%%%%%%%%%%%%%%%%%%%%%%%
v_d = u(1);   % the global voltage should be input
v_q = u(2);
P_ref = u(3);
Q_ref = u(4);

Idc=1;
C_bus=4.8e-3;

f   = 60;
w_g = 2*pi*f;

Sbase = 2e6;
Vbase = 690;
Ibase = Sbase/(sqrt(3)*Vbase);
wbase = w_g;
Zbase = Vbase^2/Sbase;
Lbase = Zbase/w_g;

% controller
Kpu=6;
Kiu=153.33;
Kpi=4.2;
Kii=0.323;
%w_f=2*pi*5;%一阶滤波器参数
%m_p=(2*pi*0.5)/Sbase;
eta=56;
alpha=325;

%滤波部分 
Lg=0.05;
Rg=7e-3;


%%%%%%%%%%%%%%%%%%%%%%%%%% Algebra %%%%%%%%%%%%%%%%%%%%%%%%
Tr_DQdq = [cos(theta) sin(theta); -sin(theta) cos(theta)]; %全局到局部
V2_dq_real = Vbase*Tr_DQdq*[v_d;v_q];
u_gd = V2_dq_real(1);
u_gq = V2_dq_real(2);

Vc_dref = v_ref;
Vc_qref = 0;
% Voltage Loop
I1_dref = Phi_vd + Kpu*(Vc_dref - u_gd);
I1_qref = Phi_vq + Kpu*(Vc_qref - u_gq);
% Current Loop
vrefd = Phi_cd - w_g*Lg*i_gq + Kpi*(I1_dref - i_gd) + Rg*i_gd;
vrefq = Phi_cq + w_g*Lg*i_gd + Kpi*(I1_qref - i_gq) + Rg*i_gq;

Pe = sqrt(3)*(u_gd*i_gd + u_gq*i_gq);
Qe = sqrt(3)*(u_gq*i_gd - u_gd*i_gq);
Vc = sqrt(u_gd^2+u_gq^2);

%%%%%%%%%%%%%%%%%%%%%%% Differential %%%%%%%%%%%%%%%%%%%%%
Tr_dqDQ = [cos(theta) -sin(theta); sin(theta) cos(theta)]; %局部到全局
is_DQ_pu = 1/Ibase*Tr_dqDQ*[i_gd;i_gq];
i_gD = is_DQ_pu(1);
i_gQ = is_DQ_pu(2);

yy = [i_gD; i_gQ];


    dotx = [
         (vrefd+w_g*Lg*i_gq-u_gd-Rg*i_gd)/Lg;
         (vrefq-w_g*Lg*i_gd-u_gq-Rg*i_gq)/Lg;
         v_ref*(eta*(Q_ref*Sbase/Vc^2-Qe/v_ref^2)+eta*alpha*(Vc^2-v_ref^2)/Vc^2);
         eta*(P_ref*Sbase/Vc^2-Pe/v_ref^2);
         Kiu*(Vc_dref-u_gd);
         Kiu*(Vc_qref-u_gq);
         Kii*(I1_dref-i_gd);
         Kii*(I1_qref-i_gq);
    ];
end


