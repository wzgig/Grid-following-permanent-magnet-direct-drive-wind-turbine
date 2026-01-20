 clc
 clear
tic

% 逆变器dq阻抗扫频相关运行程序


f_begin=1;%1
f_end=10000;
Point=20;
hp0= floor(logspace(log10(f_begin),log10(f_end),Point)) %扫描点频率
hpn=length(hp0)/5; %循环次数  20个点，每个循环注入5个频率的扰动，循环4次


 for n0=1:1:hpn
    %% 设置扰动幅值 

    funp=30;%？？
    funn=0;
    Ap1=0.01*funp;Ap2=0.03*funp;Ap3=0.05*funp;Ap4=0.07*funp;Ap5=0.09*funp;%幅值越大越精确，但过大会影响稳态工作点
    An1=0.01*funn;An2=0.03*funn;An3=0.05*funn;An4=0.07*funn;An5=0.09*funn;
    hp1=hp0(n0);hp2=hp0(n0+hpn);hp3=hp0(n0+2*hpn);hp4=hp0(n0+3*hpn);hp5=hp0(n0+4*hpn);
   fc0=[hp1,hp2,hp3,hp4,hp5]

     %% 仿真运行
    sim(['GFL_inverter.slx']);

  %% 阻抗计算
   for x=1:1:length(fc0)

         n=fc0(x)+1;

    FFTread_Vp;%算出一次
   end
 
    funp=0;
    funn=30;
    Ap1=0.01*funp;Ap2=0.03*funp;Ap3=0.05*funp;Ap4=0.07*funp;Ap5=0.09*funp;%幅值越大越精确，但过大会影响稳态工作点
    An1=0.01*funn;An2=0.03*funn;An3=0.05*funn;An4=0.07*funn;An5=0.09*funn;
    hp1=hp0(n0);hp2=hp0(n0+hpn);hp3=hp0(n0+2*hpn);hp4=hp0(n0+3*hpn);hp5=hp0(n0+4*hpn);
    fc0=[hp1,hp2,hp3,hp4,hp5]

     sim(['GFL_inverter.slx']);

 for x=1:1:length(fc0)

    n=fc0(x)+1;
    FFTread_Vn;%算出一次
  
 end
    %% 阻抗计算
len=length(fc0);
Zdq = zeros(2, 2, len); 
for x=1:1:length(fc0)
    vdqj = [vdjp(x) vdjn(x); vqjp(x) vqjn(x)];
    idqj = [idjp(x) idjn(x); iqjp(x) iqjn(x)];
    Zdq(:, :, x) = -vdqj / idqj; % 计算并存储 Z 矩阵
    GM_Zdd_scan_1(x)=GM(Zdq(1,1, x));
    PM_Zdd_scan_1(x)=PM(Zdq(1,1, x));

    GM_Zdq_scan_1(x)=GM(Zdq(1,2, x));
    PM_Zdq_scan_1(x)=PM(Zdq(1,2, x));

    GM_Zqd_scan_1(x)=GM(Zdq(2,1, x));
    PM_Zqd_scan_1(x)=PM(Zdq(2,1, x));

    GM_Zqq_scan_1(x)=GM(Zdq(2,2, x));
    PM_Zqq_scan_1(x)=PM(Zdq(2,2, x));


    
end   
     
   
% % 四个窗口
    figure(1)
    subplot(2,1,1)
    semilogx(fc0,GM_Zdd_scan_1,'r+','linewidth',1);hold on;
    subplot(2,1,2)
    semilogx(fc0,PM_Zdd_scan_1,'r+','linewidth',1);hold on;


    figure(2)
    subplot(2,1,1)
    semilogx(fc0,GM_Zdq_scan_1,'r+','linewidth',1);hold on;
    subplot(2,1,2)
    semilogx(fc0,PM_Zdq_scan_1,'r+','linewidth',1);hold on;
    
    figure(3)
    subplot(2,1,1)
    semilogx(fc0,GM_Zqd_scan_1,'r+','linewidth',1);hold on;
    subplot(2,1,2)
    semilogx(fc0,PM_Zqd_scan_1,'r+','linewidth',1);hold on;

    figure(4)
    subplot(2,1,1)
    semilogx(fc0,GM_Zqq_scan_1,'r+','linewidth',1);hold on;
    subplot(2,1,2)
    semilogx(fc0,PM_Zqq_scan_1,'r+','linewidth',1);hold on;


end 

fprintf('Run Time（Second）：');
toc
fprintf('***************END***************');